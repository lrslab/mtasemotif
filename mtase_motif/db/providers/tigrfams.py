from __future__ import annotations

import gzip
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from mtase_motif.config import DbLayout
from mtase_motif.db.base import Provider
from mtase_motif.db.manifest import ProviderArtifact, ProviderManifest
from mtase_motif.util import MtaseMotifError, run_cmd, sha256_file, which


@dataclass(frozen=True)
class ModelSelectors:
    exact: set[str]
    patterns: list[re.Pattern[str]]


class TigrfamsProvider(Provider):
    name = "tigrfams"

    def __init__(self, layout: DbLayout) -> None:
        super().__init__(layout)
        self._dir = layout.tigrfams_dir

    def fetch(self, *, source: Optional[Path] = None) -> ProviderManifest:
        if source is None:
            raise MtaseMotifError(
                "tigrfams fetch is not implemented for remote downloads yet.\n"
                "Use --source to point at a local TIGRFAMs HMM install."
            )
        source = source.expanduser().resolve()
        if not source.exists():
            raise MtaseMotifError(f"tigrfams --source does not exist: {source}")

        self._dir.mkdir(parents=True, exist_ok=True)
        src_file = _resolve_tigrfams_source(source)
        dest = self._dir / ("TIGRFAMs.hmm.gz" if src_file.name.endswith(".gz") else "TIGRFAMs.hmm")
        if dest.exists():
            dest.unlink()
        shutil.copy2(src_file, dest)

        return ProviderManifest(
            provider=self.name,
            version="unknown",
            license="TIGRFAMs (NCBI)",
            artifacts=[
                ProviderArtifact(path=str(dest.relative_to(self.layout.db_dir)), sha256=sha256_file(dest))
            ],
            notes="Raw TIGRFAMs file copied; subset build + hmmpress pending.",
        )

    def index(self) -> ProviderManifest:
        tigrfams_file = _find_tigrfams_file(self._dir)
        selectors = _load_model_selectors()
        subset_path = self._dir / "subset.hmm"
        included = _write_subset_hmm(tigrfams_file, subset_path, selectors)

        hmmpress = which("hmmpress")
        if hmmpress is None:
            raise MtaseMotifError("hmmpress not found in PATH (required for tigrfams indexing)")
        run_cmd([hmmpress, str(subset_path)], check=True, capture_output=False)

        artifacts: list[ProviderArtifact] = [
            ProviderArtifact(
                path=str(tigrfams_file.relative_to(self.layout.db_dir)),
                sha256=sha256_file(tigrfams_file),
            ),
            ProviderArtifact(
                path=str(subset_path.relative_to(self.layout.db_dir)),
                sha256=sha256_file(subset_path),
            ),
        ]
        for ext in (".h3f", ".h3i", ".h3m", ".h3p"):
            idx = Path(str(subset_path) + ext)
            if idx.exists():
                artifacts.append(
                    ProviderArtifact(
                        path=str(idx.relative_to(self.layout.db_dir)),
                        sha256=sha256_file(idx),
                    )
                )

        return ProviderManifest(
            provider=self.name,
            version="unknown",
            license="TIGRFAMs (NCBI)",
            artifacts=artifacts,
            notes=f"Built subset.hmm ({included} models) and ran hmmpress.",
        )

    def status(self) -> dict[str, object]:
        return {"provider": self.name, "dir": str(self._dir), "present": self._dir.exists()}

    def validate(self) -> None:
        if not self._dir.exists():
            raise MtaseMotifError("tigrfams not fetched")
        _find_tigrfams_file(self._dir)


def _resolve_tigrfams_source(source: Path) -> Path:
    if source.is_file():
        return source
    if source.is_dir():
        for cand in (
            source / "TIGRFAMs.hmm.gz",
            source / "TIGRFAMs.hmm",
            source / "TIGRFAMs_15.0_HMM.LIB.gz",
            source / "TIGRFAMs_15.0_HMM.LIB",
        ):
            if cand.exists():
                return cand
        for pattern in ("TIGRFAMs_*HMM.LIB.gz", "TIGRFAMs_*HMM.LIB", "*.HMM.LIB.gz", "*.HMM.LIB"):
            matches = sorted(source.glob(pattern))
            if matches:
                return matches[-1]
        raise MtaseMotifError(f"tigrfams --source dir missing TIGRFAMs HMM file: {source}")
    raise MtaseMotifError(f"tigrfams --source must be a file or dir: {source}")


def _find_tigrfams_file(tigr_dir: Path) -> Path:
    for cand in (
        tigr_dir / "TIGRFAMs.hmm.gz",
        tigr_dir / "TIGRFAMs.hmm",
        tigr_dir / "TIGRFAMs_15.0_HMM.LIB.gz",
        tigr_dir / "TIGRFAMs_15.0_HMM.LIB",
    ):
        if cand.exists():
            return cand
    matches = sorted(tigr_dir.glob("*.HMM.LIB*"))
    if matches:
        return matches[-1]
    raise MtaseMotifError(f"Missing TIGRFAMs HMM file in: {tigr_dir}")


def _load_model_selectors() -> ModelSelectors:
    models_path = Path(__file__).resolve().parents[1] / "data" / "tigrfams_models.txt"
    if not models_path.exists():
        raise MtaseMotifError(f"Missing tigrfams_models.txt: {models_path}")

    exact: set[str] = set()
    patterns: list[re.Pattern[str]] = []
    for raw in models_path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("re:"):
            expr = line[3:].strip()
            if not expr:
                continue
            patterns.append(re.compile(expr, flags=re.IGNORECASE))
        else:
            exact.add(line)

    if not exact and not patterns:
        raise MtaseMotifError(f"tigrfams model selector file is empty: {models_path}")
    return ModelSelectors(exact=exact, patterns=patterns)


def _open_text_maybe_gz(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def _write_subset_hmm(tigrfams_hmm: Path, subset_hmm: Path, selectors: ModelSelectors) -> int:
    included = 0
    record: list[str] = []
    name: Optional[str] = None
    acc: Optional[str] = None
    desc: Optional[str] = None

    subset_hmm.parent.mkdir(parents=True, exist_ok=True)
    with _open_text_maybe_gz(tigrfams_hmm) as fin, subset_hmm.open("w") as fout:
        for line in fin:
            record.append(line)
            if line.startswith("NAME  "):
                name = line.split(None, 1)[1].strip()
            elif line.startswith("ACC   "):
                acc = line.split(None, 1)[1].strip()
            elif line.startswith("DESC  "):
                desc = line.split(None, 1)[1].strip()
            elif line.strip() == "//":
                if _matches_selectors(name=name, acc=acc, desc=desc, selectors=selectors):
                    fout.writelines(record)
                    included += 1
                record = []
                name = None
                acc = None
                desc = None

    if included == 0:
        raise MtaseMotifError(
            f"No TIGRFAM models matched tigrfams_models.txt; wrote empty subset: {subset_hmm}"
        )
    return included


def _matches_selectors(
    *,
    name: Optional[str],
    acc: Optional[str],
    desc: Optional[str],
    selectors: ModelSelectors,
) -> bool:
    fields: list[str] = []
    if name:
        fields.append(name)
    if acc:
        fields.append(acc)
        fields.append(acc.split(".", 1)[0])

    if selectors.exact and any(field in selectors.exact for field in fields):
        return True

    text = " ".join([name or "", acc or "", desc or ""])
    for pattern in selectors.patterns:
        if pattern.search(text):
            return True
    return False
