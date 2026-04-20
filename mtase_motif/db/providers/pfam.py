from __future__ import annotations

import gzip
import shutil
import urllib.request
from pathlib import Path
import re
from typing import Optional, Sequence

from mtase_motif.config import DbLayout
from mtase_motif.db.base import Provider
from mtase_motif.db.manifest import ProviderArtifact, ProviderManifest
from mtase_motif.util import MtaseMotifError, run_cmd, sha256_file, which


REMOTE_PFAM_HMM_URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}?annotation=hmm"


class PfamProvider(Provider):
    name = "pfam"

    def __init__(self, layout: DbLayout) -> None:
        super().__init__(layout)
        self._dir = layout.pfam_dir

    def fetch(self, *, source: Optional[Path] = None) -> ProviderManifest:
        self._dir.mkdir(parents=True, exist_ok=True)
        if source is None:
            models = _load_model_list()
            subset_path, urls, fetched_models = _download_remote_subset_hmm(self._dir, models)
            return ProviderManifest(
                provider=self.name,
                version=_describe_subset_versions(subset_path),
                license="Pfam",
                urls=urls,
                artifacts=[
                    ProviderArtifact(
                        path=str(subset_path.relative_to(self.layout.db_dir)),
                        sha256=sha256_file(subset_path),
                    )
                ],
                notes=(
                    f"Downloaded curated Pfam subset ({len(fetched_models)} models) from the "
                    "official InterPro API; hmmpress pending."
                ),
            )

        source = source.expanduser().resolve()
        if not source.exists():
            raise MtaseMotifError(f"pfam --source does not exist: {source}")

        src_file = _resolve_pfam_source(source)
        dest = self._dir / ("Pfam-A.hmm.gz" if src_file.name.endswith(".gz") else "Pfam-A.hmm")
        if dest.exists():
            dest.unlink()
        shutil.copy2(src_file, dest)

        sha256 = sha256_file(dest)
        relpath = str(dest.relative_to(self.layout.db_dir))
        return ProviderManifest(
            provider=self.name,
            version="unknown",
            license="Pfam",
            artifacts=[ProviderArtifact(path=relpath, sha256=sha256)],
            notes="Raw Pfam file copied; subset build + hmmpress pending.",
        )

    def index(self) -> ProviderManifest:
        subset_path = self._dir / "subset.hmm"
        pfam_file = _find_pfam_a_file(self._dir)
        if pfam_file is not None:
            models = _load_model_list()
            included = _write_subset_hmm(pfam_file, subset_path, models)
            notes = f"Built subset.hmm ({included} models) from {pfam_file.name} and ran hmmpress."
        elif subset_path.exists():
            included = _count_hmm_records(subset_path)
            notes = f"Ran hmmpress on existing subset.hmm ({included} models)."
        else:
            raise MtaseMotifError(
                f"Missing Pfam inputs in: {self._dir}\n"
                "Run: mtase-motif db fetch pfam\n"
                "(or use --source <Pfam-A.hmm(.gz)> for offline mode)"
            )

        hmmpress = which("hmmpress")
        if hmmpress is None:
            raise MtaseMotifError("hmmpress not found in PATH (required for pfam indexing)")
        run_cmd([hmmpress, "-f", str(subset_path)], check=True, capture_output=False)

        artifacts: list[ProviderArtifact] = []
        if pfam_file is not None:
            artifacts.append(
                ProviderArtifact(
                    path=str(pfam_file.relative_to(self.layout.db_dir)),
                    sha256=sha256_file(pfam_file),
                )
            )
        artifacts.append(
            ProviderArtifact(
                path=str(subset_path.relative_to(self.layout.db_dir)),
                sha256=sha256_file(subset_path),
            )
        )
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
            version=_describe_subset_versions(subset_path),
            license="Pfam",
            artifacts=artifacts,
            notes=notes,
        )

    def status(self) -> dict[str, object]:
        return {"provider": self.name, "dir": str(self._dir), "present": self._dir.exists()}

    def validate(self) -> None:
        if not self._dir.exists():
            raise MtaseMotifError("pfam not fetched")


def _resolve_pfam_source(source: Path) -> Path:
    if source.is_file():
        return source
    if source.is_dir():
        for cand in (source / "Pfam-A.hmm.gz", source / "Pfam-A.hmm"):
            if cand.exists():
                return cand
        raise MtaseMotifError(f"pfam --source dir missing Pfam-A.hmm(.gz): {source}")
    raise MtaseMotifError(f"pfam --source must be a file or dir: {source}")


def _find_pfam_a_file(pfam_dir: Path) -> Optional[Path]:
    for cand in (pfam_dir / "Pfam-A.hmm.gz", pfam_dir / "Pfam-A.hmm"):
        if cand.exists():
            return cand
    return None


def _load_model_list() -> list[str]:
    models_path = Path(__file__).resolve().parents[1] / "data" / "pfam_models.txt"
    if not models_path.exists():
        raise MtaseMotifError(f"Missing pfam_models.txt: {models_path}")
    models: list[str] = []
    seen: set[str] = set()
    for raw in models_path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line in seen:
            continue
        models.append(line)
        seen.add(line)
    return models


def _open_text_maybe_gz(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def _download_remote_subset_hmm(
    pfam_dir: Path,
    models: Sequence[str],
) -> tuple[Path, list[str], list[str]]:
    urls: list[str] = []
    fetched_models: list[str] = []
    subset_hmm = pfam_dir / "subset.hmm"

    with subset_hmm.open("w") as fout:
        for model in models:
            accession = _normalize_remote_model_accession(model)
            url = REMOTE_PFAM_HMM_URL.format(accession=accession)
            record = _download_remote_hmm_record(url)
            fout.write(record)
            urls.append(url)
            fetched_models.append(accession)
    return subset_hmm, urls, fetched_models


def _normalize_remote_model_accession(model: str) -> str:
    if re.fullmatch(r"PF\d{5}(?:\.\d+)?", model):
        return model.split(".", 1)[0]
    raise MtaseMotifError(
        f"Remote Pfam fetch only supports Pfam accession IDs, got {model!r}. "
        "Use --source with a local Pfam-A.hmm(.gz) file for NAME-based subsets."
    )


def _download_remote_hmm_record(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "mtase-motif/0"})
    with urllib.request.urlopen(req) as response:
        payload = response.read()

    try:
        payload = gzip.decompress(payload)
    except OSError:
        pass

    text = payload.decode("utf-8")
    if "HMMER3/" not in text or "\n//" not in text:
        raise MtaseMotifError(f"Unexpected Pfam HMM response from {url}")
    if not text.endswith("\n"):
        text += "\n"
    return text


def _count_hmm_records(path: Path) -> int:
    with path.open() as handle:
        return sum(1 for line in handle if line.strip() == "//")


def _describe_subset_versions(path: Path) -> str:
    versions: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("ACC   "):
                versions.append(line.split(None, 1)[1].strip())
    if not versions:
        return "unknown"
    return ",".join(versions)


def _write_subset_hmm(pfam_hmm: Path, subset_hmm: Path, models: Sequence[str]) -> int:
    included = 0
    record: list[str] = []
    name: Optional[str] = None
    acc: Optional[str] = None
    model_set = set(models)

    subset_hmm.parent.mkdir(parents=True, exist_ok=True)
    with _open_text_maybe_gz(pfam_hmm) as fin, subset_hmm.open("w") as fout:
        for line in fin:
            record.append(line)
            if line.startswith("NAME  "):
                name = line.split(None, 1)[1].strip()
            elif line.startswith("ACC   "):
                acc = line.split(None, 1)[1].strip()
            elif line.strip() == "//":
                acc_base = acc.split(".", 1)[0] if acc else None
                if (name and name in model_set) or (acc and acc in model_set) or (
                    acc_base and acc_base in model_set
                ):
                    fout.writelines(record)
                    included += 1
                record = []
                name = None
                acc = None

    if included == 0:
        raise MtaseMotifError(
            f"No Pfam models matched pfam_models.txt; wrote empty subset: {subset_hmm}"
        )
    return included
