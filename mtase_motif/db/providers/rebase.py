from __future__ import annotations

import hashlib
import re
import shutil
import tempfile
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator, Optional

from mtase_motif.config import DbLayout
from mtase_motif.db.base import Provider
from mtase_motif.db.manifest import ProviderArtifact, ProviderManifest, write_provider_manifest
from mtase_motif.methylation import normalize_methylation
from mtase_motif.util import MtaseMotifError, run_cmd, sha256_file, which


REMOTE_EMBOSS_URLS: dict[str, str] = {
    "emboss_e.txt": "https://rebase.neb.com/rebase/emboss_e.txt",
    "emboss_r.txt": "https://rebase.neb.com/rebase/emboss_r.txt",
    "emboss_s.txt": "https://rebase.neb.com/rebase/emboss_s.txt",
}

REMOTE_PROTEIN_URLS: tuple[str, ...] = (
    "https://rebase.neb.com/rebase/All_REBASE_Gold_Standards_Protein.txt",
    "ftp://ftp.neb.com/pub/rebase/All_REBASE_Gold_Standards_Protein.txt",
    "ftp://ftp.neb.com/pub/rebase/All_Type_I_M_subunit_genes_Protein.txt",
    "ftp://ftp.neb.com/pub/rebase/All_Type_III_M_subunit_genes_Protein.txt",
    "ftp://ftp.neb.com/pub/rebase/All_Type_IIG_enzyme_genes_Protein.txt",
    "ftp://ftp.neb.com/pub/rebase/Type_II_methyltransferase_genes_Protein.txt",
    "ftp://ftp.neb.com/pub/rebase/Other_methyltransferase_genes_Protein.txt",
    "ftp://ftp.neb.com/pub/rebase/All_cytosine-5_methyltransferase_genes_Protein.txt",
    "ftp://ftp.neb.com/pub/rebase/All_amino-methyltransferase_genes_subtype_alpha_Protein.txt",
    "ftp://ftp.neb.com/pub/rebase/All_amino-methyltransferase_genes_subtype_beta_Protein.txt",
    "ftp://ftp.neb.com/pub/rebase/All_amino-methyltransferase_genes_subtype_gamma_Protein.txt",
    "ftp://ftp.neb.com/pub/rebase/protein_gold_seqs.txt",
    "ftp://ftp.neb.com/pub/rebase/protein_seqs.txt",
)


class RebaseProvider(Provider):
    name = "rebase"

    def __init__(self, layout: DbLayout) -> None:
        super().__init__(layout)
        self._dir = layout.rebase_dir

    def fetch(self, *, source: Optional[Path] = None) -> ProviderManifest:
        if source is None:
            with tempfile.TemporaryDirectory(prefix="_rebase_src_", dir=str(self.layout.db_dir)) as tmp:
                tmp_source = Path(tmp)
                urls, remote_note = _download_remote_rebase_inputs(tmp_source)
                return self._fetch_from_source(tmp_source, source_urls=urls, extra_note=remote_note)

        source = source.expanduser().resolve()
        if not source.exists():
            raise MtaseMotifError(f"rebase --source does not exist: {source}")
        return self._fetch_from_source(source, source_urls=None, extra_note=None)

    def _fetch_from_source(
        self,
        source: Path,
        *,
        source_urls: Optional[list[str]],
        extra_note: Optional[str],
    ) -> ProviderManifest:
        emboss_e, emboss_r, emboss_s, filename_version = _pick_emboss_files(source)
        release_version = _extract_rebase_version(emboss_e) or _extract_rebase_version(emboss_r)
        release_version = release_version or filename_version
        release_stamp = _extract_rebase_release_stamp(emboss_r)
        stamp = release_stamp or release_version or datetime.now(timezone.utc).strftime("%Y%m%d")

        dest_dir = self._dir / stamp
        raw_dir = dest_dir / "raw"
        raw_dir.mkdir(parents=True, exist_ok=True)

        e_dest = raw_dir / emboss_e.name
        r_dest = raw_dir / emboss_r.name
        s_dest = raw_dir / emboss_s.name
        shutil.copy2(emboss_e, e_dest)
        shutil.copy2(emboss_r, r_dest)
        shutil.copy2(emboss_s, s_dest)

        enzymes_tsv = dest_dir / "rebase_enzymes.tsv"
        records = _parse_emboss(e_dest, r_dest)
        _write_rebase_enzymes_tsv(enzymes_tsv, records)

        artifacts: list[ProviderArtifact] = [
            ProviderArtifact(
                path=str(e_dest.relative_to(self.layout.db_dir)),
                sha256=sha256_file(e_dest),
                url=_source_url_for_path(emboss_e, source_urls=source_urls),
            ),
            ProviderArtifact(
                path=str(r_dest.relative_to(self.layout.db_dir)),
                sha256=sha256_file(r_dest),
                url=_source_url_for_path(emboss_r, source_urls=source_urls),
            ),
            ProviderArtifact(
                path=str(s_dest.relative_to(self.layout.db_dir)),
                sha256=sha256_file(s_dest),
                url=_source_url_for_path(emboss_s, source_urls=source_urls),
            ),
            ProviderArtifact(
                path=str(enzymes_tsv.relative_to(self.layout.db_dir)),
                sha256=sha256_file(enzymes_tsv),
            ),
        ]

        proteins_dest, protein_source_urls, proteins_note = _stage_optional_proteins(source, dest_dir)
        protein_map_dest = None
        motif_proteins_dest = None
        if proteins_dest is not None:
            artifacts.append(
                ProviderArtifact(
                    path=str(proteins_dest.relative_to(self.layout.db_dir)),
                    sha256=sha256_file(proteins_dest),
                )
            )
            mapped_rows = _build_rebase_protein_map(
                records=records,
                proteins_faa=proteins_dest,
            )
            protein_map_dest = dest_dir / "rebase_protein_map.tsv"
            map_count = _write_rebase_protein_map(
                protein_map_dest,
                rows=mapped_rows,
            )
            artifacts.append(
                ProviderArtifact(
                    path=str(protein_map_dest.relative_to(self.layout.db_dir)),
                    sha256=sha256_file(protein_map_dest),
                )
            )
            motif_proteins_dest = dest_dir / "rebase_motif_proteins.faa"
            motif_protein_count = _write_rebase_motif_proteins(
                motif_proteins_dest,
                proteins_faa=proteins_dest,
                mapped_ids={row["protein_id"] for row in mapped_rows},
            )
            artifacts.append(
                ProviderArtifact(
                    path=str(motif_proteins_dest.relative_to(self.layout.db_dir)),
                    sha256=sha256_file(motif_proteins_dest),
                )
            )

        provenance_urls = list(source_urls or [])
        if not provenance_urls:
            provenance_urls.extend(
                [
                    emboss_e.resolve().as_uri(),
                    emboss_r.resolve().as_uri(),
                    emboss_s.resolve().as_uri(),
                ]
            )
            provenance_urls.extend(protein_source_urls)

        notes: list[str] = ["Parsed EMBOSS files into rebase_enzymes.tsv."]
        if source_urls:
            notes.append("Downloaded REBASE files from remote URLs.")
        if release_version:
            notes.append(f"REBASE version {release_version}.")
        if proteins_note:
            notes.append(proteins_note)
        if protein_map_dest is not None:
            notes.append(f"Built rebase_protein_map.tsv ({map_count} mapped proteins).")
        if motif_proteins_dest is not None:
            notes.append(f"Built rebase_motif_proteins.faa ({motif_protein_count} sequences).")
        if extra_note:
            notes.append(extra_note)
        if proteins_dest is None:
            notes.append(
                "No proteins FASTA found; provide rebase_proteins.faa or one or more "
                "REBASE *_Protein.txt files in --source before db index."
            )

        pm = ProviderManifest(
            provider=self.name,
            version=stamp,
            license="REBASE (NEB)",
            urls=provenance_urls,
            artifacts=artifacts,
            notes=" ".join(notes),
        )
        write_provider_manifest(dest_dir / "manifest.json", pm)
        return pm

    def index(self) -> ProviderManifest:
        dest_dir = _latest_version_dir(self._dir)
        proteins = dest_dir / "rebase_proteins.faa"
        if not proteins.exists():
            raise MtaseMotifError(
                f"Missing {proteins}. Re-run `mtase-motif db fetch rebase` "
                "or provide rebase_proteins.faa or REBASE *_Protein.txt protein dumps "
                "in --source."
            )
        mmseqs = which("mmseqs")
        makeblastdb = which("makeblastdb")

        built_index = ""
        mmseqs_dir = dest_dir / "mmseqs"
        blast_dir = dest_dir / "blast"
        if mmseqs is not None:
            mmseqs_dir.mkdir(parents=True, exist_ok=True)
            db_prefix = mmseqs_dir / "rebase_proteins"
            tmp_dir = mmseqs_dir / "tmp"
            tmp_dir.mkdir(parents=True, exist_ok=True)

            run_cmd([mmseqs, "createdb", str(proteins), str(db_prefix)], check=True, capture_output=False)
            run_cmd(
                [mmseqs, "createindex", str(db_prefix), str(tmp_dir)],
                check=True,
                capture_output=False,
            )
            built_index = "mmseqs"
        elif makeblastdb is not None:
            blast_dir.mkdir(parents=True, exist_ok=True)
            db_prefix = blast_dir / "rebase_proteins"
            run_cmd(
                [
                    makeblastdb,
                    "-in",
                    str(proteins),
                    "-dbtype",
                    "prot",
                    "-out",
                    str(db_prefix),
                ],
                check=True,
                capture_output=False,
            )
            built_index = "blast"
        else:
            raise MtaseMotifError("Neither mmseqs nor makeblastdb found in PATH (need one for indexing)")

        artifacts: list[ProviderArtifact] = []
        raw_dir = dest_dir / "raw"
        if raw_dir.exists():
            for p in sorted(raw_dir.glob("emboss_*")):
                if p.is_file():
                    artifacts.append(
                        ProviderArtifact(
                            path=str(p.relative_to(self.layout.db_dir)),
                            sha256=sha256_file(p),
                        )
                    )
        enzymes_tsv = dest_dir / "rebase_enzymes.tsv"
        protein_map_tsv = dest_dir / "rebase_protein_map.tsv"
        motif_proteins_faa = dest_dir / "rebase_motif_proteins.faa"
        cluster_labels_tsv = dest_dir / "rebase_cluster_labels.tsv"
        cluster_label_count = 0
        if enzymes_tsv.exists():
            records = _load_rebase_enzymes_tsv(enzymes_tsv)
            _write_rebase_enzymes_tsv(enzymes_tsv, records)
            artifacts.append(
                ProviderArtifact(
                    path=str(enzymes_tsv.relative_to(self.layout.db_dir)),
                    sha256=sha256_file(enzymes_tsv),
                )
            )
            mapped_rows = _build_rebase_protein_map(
                records=records,
                proteins_faa=proteins,
            )
            _write_rebase_protein_map(
                protein_map_tsv,
                rows=mapped_rows,
            )
            artifacts.append(
                ProviderArtifact(
                    path=str(protein_map_tsv.relative_to(self.layout.db_dir)),
                    sha256=sha256_file(protein_map_tsv),
                )
            )
            _write_rebase_motif_proteins(
                motif_proteins_faa,
                proteins_faa=proteins,
                mapped_ids={row["protein_id"] for row in mapped_rows},
            )
            artifacts.append(
                ProviderArtifact(
                    path=str(motif_proteins_faa.relative_to(self.layout.db_dir)),
                    sha256=sha256_file(motif_proteins_faa),
                )
            )
            cluster_label_count = _write_rebase_cluster_labels(
                cluster_labels_tsv,
                proteins_faa=proteins,
                motif_proteins_faa=motif_proteins_faa,
                protein_map_rows=mapped_rows,
                blastp_exe=which("blastp"),
            )
            artifacts.append(
                ProviderArtifact(
                    path=str(cluster_labels_tsv.relative_to(self.layout.db_dir)),
                    sha256=sha256_file(cluster_labels_tsv),
                )
            )
        artifacts.append(
            ProviderArtifact(
                path=str(proteins.relative_to(self.layout.db_dir)),
                sha256=sha256_file(proteins),
            )
        )
        if built_index == "mmseqs":
            for p in sorted(mmseqs_dir.glob("rebase_proteins*")):
                if p.is_file():
                    artifacts.append(
                        ProviderArtifact(
                            path=str(p.relative_to(self.layout.db_dir)),
                            sha256=sha256_file(p),
                        )
                    )
        elif built_index == "blast":
            for p in sorted(blast_dir.glob("rebase_proteins*")):
                if p.is_file():
                    artifacts.append(
                        ProviderArtifact(
                            path=str(p.relative_to(self.layout.db_dir)),
                            sha256=sha256_file(p),
                        )
                    )

        pm = ProviderManifest(
            provider=self.name,
            version=dest_dir.name,
            license="REBASE (NEB)",
            artifacts=artifacts,
            notes=(
                f"Built {built_index} protein index for rebase_proteins.faa and refreshed "
                f"rebase_protein_map.tsv and rebase_cluster_labels.tsv ({cluster_label_count} labeled proteins)."
            ),
        )
        write_provider_manifest(dest_dir / "manifest.json", pm)
        return pm

    def status(self) -> dict[str, object]:
        latest = ""
        if self._dir.exists():
            dirs = sorted([p for p in self._dir.iterdir() if p.is_dir()])
            if dirs:
                latest = dirs[-1].name
        return {
            "provider": self.name,
            "dir": str(self._dir),
            "present": self._dir.exists(),
            "latest_version": latest,
        }

    def validate(self) -> None:
        if not self._dir.exists():
            raise MtaseMotifError("rebase not fetched")
        latest = _latest_version_dir(self._dir)
        enzymes = latest / "rebase_enzymes.tsv"
        if not enzymes.exists():
            raise MtaseMotifError(f"Missing parsed enzymes table: {enzymes}")


def _download_remote_rebase_inputs(dest_dir: Path) -> tuple[list[str], Optional[str]]:
    urls: list[str] = []
    for filename, url in REMOTE_EMBOSS_URLS.items():
        _download_url(url, dest_dir / filename)
        urls.append(url)

    protein_urls, proteins_note = _download_remote_proteins(dest_dir)
    if protein_urls:
        urls.extend(protein_urls)
        return urls, proteins_note
    return (
        urls,
        "Could not auto-download REBASE proteins FASTA from known endpoints.",
    )


def _download_remote_proteins(dest_dir: Path) -> tuple[list[str], Optional[str]]:
    staging_dir = dest_dir / "_rebase_proteins_remote"
    out_faa = dest_dir / "rebase_proteins.faa"
    staging_dir.mkdir(parents=True, exist_ok=True)

    normalized_parts: list[Path] = []
    used_urls: list[str] = []
    for idx, url in enumerate(REMOTE_PROTEIN_URLS):
        name = Path(urllib.parse.urlparse(url).path).name or f"part_{idx}"
        raw_path = staging_dir / f"{idx:02d}_{name}"
        norm_path = staging_dir / f"{idx:02d}_{Path(name).stem}.faa"
        try:
            _download_url(url, raw_path)
            seq_count = _normalize_rebase_proteins(raw_path, norm_path)
            if seq_count > 0:
                normalized_parts.append(norm_path)
                used_urls.append(url)
        except MtaseMotifError:
            continue

    if not normalized_parts:
        shutil.rmtree(staging_dir, ignore_errors=True)
        out_faa.unlink(missing_ok=True)
        return [], None

    unique_count = _merge_normalized_fastas(normalized_parts, out_faa)
    shutil.rmtree(staging_dir, ignore_errors=True)
    return (
        used_urls,
        f"Downloaded and normalized {len(used_urls)} remote REBASE protein set(s) into "
        f"{unique_count} unique sequences.",
    )


def _download_url(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        if url.startswith("ftp://"):
            req = url
        else:
            req = urllib.request.Request(url, headers={"User-Agent": "mtase-motif/0.1"})
        with urllib.request.urlopen(req, timeout=90) as response, dest.open("wb") as out:
            shutil.copyfileobj(response, out)
    except Exception as exc:
        raise MtaseMotifError(f"Failed to download {url}: {exc}") from exc

    if not dest.exists() or dest.stat().st_size == 0:
        raise MtaseMotifError(f"Downloaded empty file from {url}")


def _source_url_for_path(path: Path, *, source_urls: Optional[list[str]]) -> Optional[str]:
    if source_urls:
        return REMOTE_EMBOSS_URLS.get(path.name)
    return path.resolve().as_uri()


def _pick_emboss_files(source: Path) -> tuple[Path, Path, Path, str]:
    if source.is_file():
        raise MtaseMotifError("rebase --source must be a directory")
    if not source.is_dir():
        raise MtaseMotifError("rebase --source must be a directory")

    role_files = {
        "e": sorted(source.glob("emboss_e.*")),
        "r": sorted(source.glob("emboss_r.*")),
        "s": sorted(source.glob("emboss_s.*")),
    }
    if not all(role_files.values()):
        raise MtaseMotifError(
            f"Missing REBASE EMBOSS files in {source}. Expected emboss_e.*, emboss_r.*, emboss_s.*"
        )

    groups: dict[str, dict[str, Path]] = {}
    for role, paths in role_files.items():
        for path in paths:
            suffix = path.name.split(".", 1)[1] if "." in path.name else ""
            groups.setdefault(suffix, {})[role] = path

    complete_groups = {
        suffix: mapping
        for suffix, mapping in groups.items()
        if {"e", "r", "s"} <= set(mapping)
    }
    if not complete_groups:
        available = ", ".join(sorted(suffix or "<no suffix>" for suffix in groups))
        raise MtaseMotifError(
            "Found REBASE EMBOSS files, but not a single coherent release with matching "
            f"`emboss_e.*`, `emboss_r.*`, and `emboss_s.*` suffixes in {source}. "
            f"Available suffix groups: {available}"
        )

    numeric_suffixes = [suffix for suffix in complete_groups if suffix.isdigit()]
    if numeric_suffixes:
        chosen_suffix = max(numeric_suffixes, key=int)
    elif len(complete_groups) == 1:
        chosen_suffix = next(iter(complete_groups))
    else:
        available = ", ".join(sorted(suffix or "<no suffix>" for suffix in complete_groups))
        raise MtaseMotifError(
            "Multiple complete non-numeric REBASE EMBOSS release groups were found in "
            f"{source}; please provide a directory containing a single release. "
            f"Complete groups: {available}"
        )

    chosen = complete_groups[chosen_suffix]
    version = chosen_suffix if chosen_suffix.isdigit() else ""
    return chosen["e"], chosen["r"], chosen["s"], version


def _extract_rebase_version(path: Path) -> str:
    for line in _head_lines(path, 120):
        m = re.search(r"REBASE version\s+([0-9]+)", line, flags=re.IGNORECASE)
        if m:
            return m.group(1)
        m = re.search(r"emboss_[ers]\.([0-9]+)", line, flags=re.IGNORECASE)
        if m:
            return m.group(1)
    return ""


def _extract_rebase_release_stamp(path: Path) -> str:
    for line in _head_lines(path, 200):
        m = re.search(r"Rich Roberts\s+([A-Za-z]{3}\s+\d{1,2}\s+\d{4})", line)
        if not m:
            continue
        raw_date = " ".join(m.group(1).split())
        try:
            dt = datetime.strptime(raw_date, "%b %d %Y")
            return dt.strftime("%Y%m%d")
        except ValueError:
            continue
    return ""


def _head_lines(path: Path, max_lines: int) -> list[str]:
    lines: list[str] = []
    with path.open(errors="ignore") as f:
        for _i, raw in enumerate(f):
            lines.append(raw.rstrip("\n"))
            if len(lines) >= max_lines:
                break
    return lines


def _stage_optional_proteins(
    source_dir: Path,
    dest_dir: Path,
) -> tuple[Optional[Path], list[str], Optional[str]]:
    proteins_sources = _find_optional_proteins_files(source_dir)
    if not proteins_sources:
        return None, [], None

    proteins_dest = dest_dir / "rebase_proteins.faa"
    source_urls = [src.resolve().as_uri() for src in proteins_sources]
    if len(proteins_sources) == 1:
        proteins_src = proteins_sources[0]
        seq_count = _normalize_rebase_proteins(proteins_src, proteins_dest)
        if seq_count <= 0:
            proteins_dest.unlink(missing_ok=True)
            raise MtaseMotifError(f"Could not parse any protein sequences from: {proteins_src}")
        note = f"Normalized proteins from {proteins_src.name} ({seq_count} sequences)."
        return proteins_dest, source_urls, note

    with tempfile.TemporaryDirectory(prefix="_rebase_proteins_", dir=str(dest_dir)) as tmp:
        staging_dir = Path(tmp)
        normalized_parts: list[Path] = []
        parsed_sources: list[Path] = []
        for idx, proteins_src in enumerate(proteins_sources):
            norm_path = staging_dir / f"{idx:02d}_{proteins_src.stem}.faa"
            seq_count = _normalize_rebase_proteins(proteins_src, norm_path)
            if seq_count <= 0:
                continue
            normalized_parts.append(norm_path)
            parsed_sources.append(proteins_src)

        if not normalized_parts:
            proteins_dest.unlink(missing_ok=True)
            raise MtaseMotifError(
                "Could not parse any protein sequences from the REBASE source files: "
                + ", ".join(str(src) for src in proteins_sources)
            )

        unique_count = _merge_normalized_fastas(normalized_parts, proteins_dest)

    note = (
        f"Merged {len(parsed_sources)} REBASE protein source files into {unique_count} unique sequences: "
        + ", ".join(src.name for src in parsed_sources)
    )
    return proteins_dest, source_urls, note


def _find_optional_proteins_files(source_dir: Path) -> list[Path]:
    if not source_dir.is_dir():
        return []

    for name in ("rebase_proteins.faa", "rebase_proteins.fasta", "rebase_proteins.fa"):
        p = source_dir / name
        if p.exists():
            return [p]

    preferred_names = (
        "All_REBASE_Gold_Standards_Protein.txt",
        "Type_II_methyltransferase_genes_Protein.txt",
        "Other_methyltransferase_genes_Protein.txt",
        "All_cytosine-5_methyltransferase_genes_Protein.txt",
        "All_amino-methyltransferase_genes_subtype_alpha_Protein.txt",
        "All_amino-methyltransferase_genes_subtype_beta_Protein.txt",
        "All_amino-methyltransferase_genes_subtype_gamma_Protein.txt",
        "protein_gold_seqs.txt",
        "protein_seqs.txt",
    )
    sources: list[Path] = []
    seen: set[Path] = set()

    for name in preferred_names:
        p = source_dir / name
        if p.exists():
            resolved = p.resolve()
            if resolved not in seen:
                sources.append(p)
                seen.add(resolved)

    for p in sorted(source_dir.glob("*_Protein.txt")):
        resolved = p.resolve()
        if resolved in seen:
            continue
        sources.append(p)
        seen.add(resolved)

    return sources


def _normalize_rebase_header_id(header: str) -> str:
    token = re.split(r"\s+", header.strip(), maxsplit=1)[0]
    if token.upper().startswith("REBASE:"):
        token = token.split(":", 1)[1]
    token = re.sub(r"[^A-Za-z0-9_.-]", "_", token)
    return token


def _normalize_rebase_proteins(source: Path, dest: Path) -> int:
    dest.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    seen_ids: set[str] = set()
    current_id: Optional[str] = None
    sequence_chunks: list[str] = []

    def flush(out_handle) -> None:
        nonlocal count, current_id, sequence_chunks
        if current_id is None:
            return

        sequence = "".join(sequence_chunks)
        rec_id = current_id
        current_id = None
        sequence_chunks = []

        if not sequence:
            return
        if rec_id in seen_ids:
            return
        seen_ids.add(rec_id)

        out_handle.write(f">{rec_id}\n")
        for i in range(0, len(sequence), 60):
            out_handle.write(sequence[i : i + 60] + "\n")
        count += 1

    with source.open(errors="ignore") as fin, dest.open("w") as fout:
        for raw in fin:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush(fout)
                current_id = _normalize_rebase_header_id(line[1:]) or f"rebase_seq_{count + 1}"
                continue
            if current_id is None:
                continue
            if line == "<>":
                continue
            aa = "".join(ch for ch in line.upper() if "A" <= ch <= "Z")
            if aa:
                sequence_chunks.append(aa)

        flush(fout)

    return count


def _merge_normalized_fastas(sources: list[Path], dest: Path) -> int:
    dest.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    seen_ids: set[str] = set()

    with dest.open("w") as fout:
        for source in sources:
            current_id: Optional[str] = None
            sequence_chunks: list[str] = []

            def flush() -> None:
                nonlocal count, current_id, sequence_chunks
                if current_id is None:
                    return
                sequence = "".join(sequence_chunks)
                rec_id = current_id
                current_id = None
                sequence_chunks = []
                if not sequence or rec_id in seen_ids:
                    return
                seen_ids.add(rec_id)
                fout.write(f">{rec_id}\n")
                for i in range(0, len(sequence), 60):
                    fout.write(sequence[i : i + 60] + "\n")
                count += 1

            with source.open(errors="ignore") as fin:
                for raw in fin:
                    line = raw.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        flush()
                        current_id = line[1:].strip()
                        continue
                    if current_id is None:
                        continue
                    sequence_chunks.append(line)
                flush()

    return count


def _parse_emboss(emboss_e: Path, emboss_r: Path) -> list[dict[str, str]]:
    enzyme_to_pattern: dict[str, str] = {}
    for raw in emboss_e.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        enzyme, pattern = parts[0], parts[1]
        enzyme_to_pattern[enzyme] = pattern.upper()

    records: list[dict[str, str]] = []
    lines = emboss_r.read_text().splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].rstrip("\n")
        if not line or line.startswith("#"):
            i += 1
            continue
        enzyme_id = line.strip()
        if i + 7 >= len(lines):
            break
        organism = lines[i + 1].strip()
        isoschizomers = lines[i + 2].strip()
        methylation = lines[i + 3].strip()
        try:
            nrefs = int(lines[i + 6].strip())
        except Exception:
            nrefs = 0
        refs_start = i + 7
        refs_end = min(refs_start + nrefs, len(lines))
        i = refs_end
        while i < len(lines) and lines[i].strip() != "//":
            i += 1
        if i < len(lines) and lines[i].strip() == "//":
            i += 1

        pattern = enzyme_to_pattern.get(enzyme_id, "")
        if isoschizomers:
            aliases = [x.strip() for x in isoschizomers.split(",") if x.strip()]
            names = ",".join([enzyme_id] + aliases)
        else:
            names = enzyme_id
        records.append(
            {
                "enzyme_id": enzyme_id,
                "names": names,
                "type": "unknown",
                "recognition_seq_iupac": pattern,
                "methylation": normalize_methylation(methylation),
                "methylation_raw": methylation,
                "evidence": "rebase_emboss",
                "organism": organism,
            }
        )

    return records


def _write_rebase_enzymes_tsv(path: Path, records: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "enzyme_id",
        "names",
        "type",
        "recognition_seq_iupac",
        "methylation",
        "methylation_raw",
        "evidence",
    ]
    with path.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for rec in records:
            f.write("\t".join(rec.get(c, "") for c in cols) + "\n")


def _load_rebase_enzymes_tsv(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    lines = path.read_text().splitlines()
    if not lines:
        return records

    header = lines[0].split("\t")
    for raw in lines[1:]:
        if not raw:
            continue
        parts = raw.split("\t")
        rec = {header[i]: parts[i] if i < len(parts) else "" for i in range(len(header))}
        methylation_raw = rec.get("methylation_raw", "") or rec.get("methylation", "")
        rec["methylation_raw"] = methylation_raw
        rec["methylation"] = normalize_methylation(rec.get("methylation", ""))
        records.append(rec)
    return records


def _write_rebase_protein_map(path: Path, *, rows: list[dict[str, str]]) -> int:
    cols = [
        "protein_id",
        "enzyme_id",
        "matched_name",
        "motif_iupac",
        "methylation",
        "methylation_raw",
        "mapping_method",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for row in rows:
            f.write("\t".join(row.get(col, "") for col in cols) + "\n")
    return len(rows)


def _build_rebase_protein_map(
    *,
    records: list[dict[str, str]],
    proteins_faa: Path,
) -> list[dict[str, str]]:
    alias_index: dict[str, list[tuple[dict[str, str], str, str]]] = {}
    for rec in records:
        names: list[str] = []
        enzyme_id = rec.get("enzyme_id", "").strip()
        if enzyme_id:
            names.append(enzyme_id)
        extra_names = rec.get("names", "")
        if extra_names:
            names.extend([part.strip() for part in extra_names.split(",") if part.strip()])

        seen_aliases: set[tuple[str, str]] = set()
        for alias in names:
            for key, method in _rebase_name_keys(alias):
                pair = (key, method)
                if pair in seen_aliases:
                    continue
                seen_aliases.add(pair)
                alias_index.setdefault(key, []).append((rec, alias, method))

    rows: list[dict[str, str]] = []
    for protein_id in _iter_rebase_protein_ids(proteins_faa):
        match = _match_rebase_protein_id(protein_id, alias_index)
        if match is None:
            continue
        rec, matched_name, method = match
        rows.append(
            {
                "protein_id": protein_id,
                "enzyme_id": rec.get("enzyme_id", ""),
                "matched_name": matched_name,
                "motif_iupac": rec.get("recognition_seq_iupac", ""),
                "methylation": normalize_methylation(rec.get("methylation", "")),
                "methylation_raw": rec.get("methylation_raw", "") or rec.get("methylation", ""),
                "mapping_method": method,
            }
        )
    rows.sort(key=lambda row: (row["protein_id"], row["enzyme_id"], row["mapping_method"]))
    return rows


def _iter_rebase_protein_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open(errors="ignore") as f:
        for raw in f:
            if not raw.startswith(">"):
                continue
            protein_id = _normalize_rebase_header_id(raw[1:])
            if protein_id:
                ids.append(protein_id)
    return ids


def _write_rebase_motif_proteins(
    path: Path,
    *,
    proteins_faa: Path,
    mapped_ids: set[str],
) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    keep = False
    with proteins_faa.open(errors="ignore") as fin, path.open("w") as fout:
        for raw in fin:
            if raw.startswith(">"):
                protein_id = _normalize_rebase_header_id(raw[1:])
                keep = protein_id in mapped_ids
                if keep:
                    fout.write(f">{protein_id}\n")
                    count += 1
                continue
            if keep:
                fout.write(raw)
    return count


def _write_rebase_cluster_labels(
    path: Path,
    *,
    proteins_faa: Path,
    motif_proteins_faa: Path,
    protein_map_rows: list[dict[str, str]],
    blastp_exe: Optional[str],
) -> int:
    direct_labels, pending_clusters = _build_rebase_duplicate_clusters(
        proteins_faa=proteins_faa,
        protein_map_rows=protein_map_rows,
    )
    neighbor_labels: dict[str, dict[str, str]] = {}
    if pending_clusters and blastp_exe is not None and motif_proteins_faa.exists():
        neighbor_labels = _infer_rebase_duplicate_neighbor_labels(
            pending_clusters=pending_clusters,
            motif_proteins_faa=motif_proteins_faa,
            protein_map_rows=protein_map_rows,
            blastp_exe=blastp_exe,
        )

    labels_by_cluster = {**direct_labels, **neighbor_labels}
    cols = [
        "protein_id",
        "cluster_id",
        "representative_id",
        "cluster_size",
        "motif_iupac",
        "methylation",
        "confidence",
        "label_method",
        "best_subject_id",
        "best_evalue",
        "best_bits",
        "best_pident",
        "best_qcov",
        "best_tcov",
        "support_count",
    ]
    count = 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for protein_id, seq in _iter_fasta_records(proteins_faa):
            cluster_id = _rebase_sequence_cluster_id(seq)
            label = labels_by_cluster.get(cluster_id)
            if label is None:
                continue
            row = {"protein_id": protein_id, **label}
            f.write("\t".join(row.get(col, "") for col in cols) + "\n")
            count += 1
    return count


def _build_rebase_duplicate_clusters(
    *,
    proteins_faa: Path,
    protein_map_rows: list[dict[str, str]],
) -> tuple[dict[str, dict[str, str]], list[dict[str, str]]]:
    mapped_by_id = {row["protein_id"]: row for row in protein_map_rows if row.get("protein_id")}
    clusters: dict[str, dict[str, object]] = {}
    for protein_id, seq in _iter_fasta_records(proteins_faa):
        cluster_id = _rebase_sequence_cluster_id(seq)
        cluster = clusters.setdefault(
            cluster_id,
            {
                "cluster_id": cluster_id,
                "representative_id": protein_id,
                "sequence": seq,
                "cluster_size": 0,
                "mapped_groups": {},
                "mapped_subjects": {},
            },
        )
        cluster["cluster_size"] = int(cluster["cluster_size"]) + 1
        mapped = mapped_by_id.get(protein_id)
        if mapped is None or not mapped.get("motif_iupac"):
            continue
        key = (mapped.get("motif_iupac", ""), normalize_methylation(mapped.get("methylation", "")))
        mapped_groups = cluster["mapped_groups"]
        mapped_groups[key] = int(mapped_groups.get(key, 0)) + 1
        mapped_subjects = cluster["mapped_subjects"]
        mapped_subjects.setdefault(key, protein_id)

    direct_labels: dict[str, dict[str, str]] = {}
    pending_clusters: list[dict[str, str]] = []
    for cluster_id, cluster in clusters.items():
        cluster_size = int(cluster["cluster_size"])
        if cluster_size < 2:
            continue

        mapped_groups = cluster["mapped_groups"]
        if len(mapped_groups) == 1:
            (motif_iupac, methylation), support_count = next(iter(mapped_groups.items()))
            if motif_iupac:
                mapped_subjects = cluster["mapped_subjects"]
                direct_labels[cluster_id] = {
                    "cluster_id": cluster_id,
                    "representative_id": str(cluster["representative_id"]),
                    "cluster_size": str(cluster_size),
                    "motif_iupac": motif_iupac,
                    "methylation": methylation,
                    "confidence": "high",
                    "label_method": "direct_duplicate_consensus",
                    "best_subject_id": str(mapped_subjects.get((motif_iupac, methylation), "")),
                    "best_evalue": "",
                    "best_bits": "",
                    "best_pident": "",
                    "best_qcov": "",
                    "best_tcov": "",
                    "support_count": str(support_count),
                }
            continue

        if cluster_size >= 20:
            pending_clusters.append(
                {
                    "cluster_id": cluster_id,
                    "representative_id": str(cluster["representative_id"]),
                    "sequence": str(cluster["sequence"]),
                    "cluster_size": str(cluster_size),
                }
            )
    return direct_labels, pending_clusters


def _infer_rebase_duplicate_neighbor_labels(
    *,
    pending_clusters: list[dict[str, str]],
    motif_proteins_faa: Path,
    protein_map_rows: list[dict[str, str]],
    blastp_exe: str,
) -> dict[str, dict[str, str]]:
    if not pending_clusters:
        return {}

    protein_map_by_id = {row["protein_id"]: row for row in protein_map_rows if row.get("protein_id")}
    hits_by_cluster = _blast_duplicate_cluster_representatives(
        pending_clusters=pending_clusters,
        motif_proteins_faa=motif_proteins_faa,
        blastp_exe=blastp_exe,
    )

    labels: dict[str, dict[str, str]] = {}
    for cluster in pending_clusters:
        cluster_id = cluster["cluster_id"]
        label = _choose_duplicate_neighbor_label(
            cluster=cluster,
            hits=hits_by_cluster.get(cluster_id, []),
            protein_map_by_id=protein_map_by_id,
        )
        if label is not None:
            labels[cluster_id] = label
    return labels


def _blast_duplicate_cluster_representatives(
    *,
    pending_clusters: list[dict[str, str]],
    motif_proteins_faa: Path,
    blastp_exe: str,
) -> dict[str, list[dict[str, object]]]:
    if not pending_clusters:
        return {}

    with tempfile.TemporaryDirectory(prefix="_rebase_cluster_blast_") as tmpdir:
        query_faa = Path(tmpdir) / "cluster_representatives.faa"
        with query_faa.open("w") as f:
            for cluster in pending_clusters:
                f.write(f">{cluster['cluster_id']}\n{cluster['sequence']}\n")
        proc = run_cmd(
            [
                blastp_exe,
                "-query",
                str(query_faa),
                "-subject",
                str(motif_proteins_faa),
                "-max_target_seqs",
                "10",
                "-max_hsps",
                "1",
                "-outfmt",
                "6 qseqid sseqid evalue bitscore pident qcovs slen length",
            ],
            check=True,
            capture_output=True,
        )

    hits_by_cluster: dict[str, list[dict[str, object]]] = {}
    for raw in proc.stdout.splitlines():
        if not raw.strip():
            continue
        parts = raw.split("\t")
        if len(parts) != 8:
            continue
        try:
            evalue = float(parts[2])
            bits = float(parts[3])
            pident = float(parts[4])
            qcov = float(parts[5])
            slen = float(parts[6])
            alen = float(parts[7])
        except ValueError:
            continue
        tcov = 0.0 if slen <= 0 else (100.0 * alen / slen)
        hits_by_cluster.setdefault(parts[0], []).append(
            {
                "subject_id": parts[1],
                "evalue": evalue,
                "bits": bits,
                "pident": pident,
                "qcov": qcov,
                "tcov": tcov,
            }
        )
    return hits_by_cluster


def _choose_duplicate_neighbor_label(
    *,
    cluster: dict[str, str],
    hits: list[dict[str, object]],
    protein_map_by_id: dict[str, dict[str, str]],
) -> Optional[dict[str, str]]:
    groups: dict[tuple[str, str], dict[str, object]] = {}
    for hit in hits:
        subject_id = str(hit["subject_id"])
        mapped = protein_map_by_id.get(subject_id)
        if mapped is None or not mapped.get("motif_iupac"):
            continue

        evalue = float(hit["evalue"])
        bits = float(hit["bits"])
        pident = float(hit["pident"])
        qcov = float(hit["qcov"])
        tcov = float(hit["tcov"])
        if evalue > 1e-10 or bits < 60.0 or qcov < 75.0 or tcov < 40.0:
            continue

        key = (
            mapped.get("motif_iupac", ""),
            normalize_methylation(mapped.get("methylation", "")),
        )
        group = groups.setdefault(
            key,
            {
                "score": 0.0,
                "support_count": 0,
                "best_subject_id": "",
                "best_evalue": float("inf"),
                "best_bits": 0.0,
                "best_pident": 0.0,
                "best_qcov": 0.0,
                "best_tcov": 0.0,
            },
        )
        group["score"] = float(group["score"]) + bits
        group["support_count"] = int(group["support_count"]) + 1
        best_tuple = (
            evalue,
            -bits,
            -pident,
            -qcov,
            -tcov,
        )
        current_tuple = (
            float(group["best_evalue"]),
            -float(group["best_bits"]),
            -float(group["best_pident"]),
            -float(group["best_qcov"]),
            -float(group["best_tcov"]),
        )
        if best_tuple < current_tuple:
            group["best_subject_id"] = subject_id
            group["best_evalue"] = evalue
            group["best_bits"] = bits
            group["best_pident"] = pident
            group["best_qcov"] = qcov
            group["best_tcov"] = tcov

    if not groups:
        return None

    ranked = sorted(
        groups.items(),
        key=lambda item: (
            -float(item[1]["score"]),
            -int(item[1]["support_count"]),
            float(item[1]["best_evalue"]),
        ),
    )
    best_key, best = ranked[0]
    second_score = float(ranked[1][1]["score"]) if len(ranked) > 1 else 0.0

    best_bits = float(best["best_bits"])
    best_pident = float(best["best_pident"])
    best_qcov = float(best["best_qcov"])
    best_tcov = float(best["best_tcov"])
    support_count = int(best["support_count"])
    accept = False
    confidence = ""
    if support_count >= 2 and float(best["score"]) >= max(140.0, second_score * 1.35):
        accept = True
        confidence = "medium" if best_bits >= 150.0 and best_qcov >= 80.0 and best_tcov >= 70.0 else "low"
    elif (
        support_count >= 1
        and best_bits >= 250.0
        and best_pident >= 50.0
        and best_qcov >= 80.0
        and best_tcov >= 80.0
        and float(best["score"]) >= max(250.0, second_score * 1.5)
    ):
        accept = True
        confidence = "high" if best_bits >= 500.0 or best_pident >= 80.0 else "medium"
    if not accept or not best_key[0]:
        return None

    return {
        "cluster_id": cluster["cluster_id"],
        "representative_id": cluster["representative_id"],
        "cluster_size": cluster["cluster_size"],
        "motif_iupac": best_key[0],
        "methylation": best_key[1],
        "confidence": confidence,
        "label_method": "neighbor_duplicate_consensus",
        "best_subject_id": str(best["best_subject_id"]),
        "best_evalue": f"{float(best['best_evalue']):g}",
        "best_bits": f"{best_bits:g}",
        "best_pident": f"{best_pident:g}",
        "best_qcov": f"{best_qcov:g}",
        "best_tcov": f"{best_tcov:g}",
        "support_count": str(support_count),
    }


def _iter_fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    current_id: Optional[str] = None
    chunks: list[str] = []

    def flush() -> Optional[tuple[str, str]]:
        nonlocal current_id, chunks
        if current_id is None:
            return None
        record = (current_id, "".join(chunks))
        current_id = None
        chunks = []
        return record

    with path.open(errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                record = flush()
                if record is not None:
                    yield record
                current_id = _normalize_rebase_header_id(line[1:])
                continue
            if current_id is not None:
                chunks.append(line.rstrip("*"))
    record = flush()
    if record is not None:
        yield record


def _rebase_sequence_cluster_id(sequence: str) -> str:
    return hashlib.sha1(sequence.encode()).hexdigest()[:16]


def _match_rebase_protein_id(
    protein_id: str,
    alias_index: dict[str, list[tuple[dict[str, str], str, str]]],
) -> Optional[tuple[dict[str, str], str, str]]:
    for key, key_method in _rebase_name_keys(protein_id):
        matches = alias_index.get(key, [])
        if not matches:
            continue

        ranked = sorted(
            matches,
            key=lambda item: (
                0 if item[0].get("recognition_seq_iupac", "") else 1,
                item[0].get("enzyme_id", ""),
                item[1],
            ),
        )
        best = ranked[0]
        best_motif = (
            best[0].get("recognition_seq_iupac", ""),
            best[0].get("methylation", ""),
        )
        conflicting = [
            item
            for item in ranked[1:]
            if (
                item[0].get("recognition_seq_iupac", ""),
                item[0].get("methylation", ""),
            )
            != best_motif
        ]
        if conflicting:
            continue
        return best[0], best[1], f"alias_{key_method}"
    return None


def _rebase_name_keys(value: str) -> list[tuple[str, str]]:
    raw = value.strip()
    if not raw:
        return []

    variants = [
        ("exact", raw),
        ("strip_prefix", re.sub(r"^[A-Za-z]\.", "", raw)),
        ("strip_suffix", raw[:-1] if raw.endswith("P") else raw),
    ]
    stripped = re.sub(r"^[A-Za-z]\.", "", raw)
    variants.append(("strip_prefix_suffix", stripped[:-1] if stripped.endswith("P") else stripped))

    keys: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for method, variant in variants:
        key = re.sub(r"[^A-Za-z0-9]", "", variant).lower()
        if not key:
            continue
        pair = (key, method)
        if pair in seen:
            continue
        seen.add(pair)
        keys.append(pair)
    return keys


def _latest_version_dir(rebase_dir: Path) -> Path:
    if not rebase_dir.exists():
        raise MtaseMotifError(f"Missing rebase dir: {rebase_dir}")
    dirs = sorted([p for p in rebase_dir.iterdir() if p.is_dir()])
    if not dirs:
        raise MtaseMotifError(f"No rebase versions found under: {rebase_dir}")
    return dirs[-1]
