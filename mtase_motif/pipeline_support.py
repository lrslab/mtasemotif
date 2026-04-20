from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Iterable, Optional

from mtase_motif.util import MtaseMotifError


STRUCTURE_EXTS = (
    ".pdb",
    ".cif",
    ".mmcif",
    ".pdb.gz",
    ".cif.gz",
    ".mmcif.gz",
)


def write_placeholder_gff(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("# No gene calls produced (proteins were provided via --proteins)\n")


def extract_candidate_proteins(*, proteins_faa: Path, candidates_tsv: Path, out_faa: Path) -> None:
    candidate_ids = load_candidate_ids(candidates_tsv)
    out_faa.parent.mkdir(parents=True, exist_ok=True)

    with out_faa.open("w") as out:
        for seq_id, seq in iter_fasta(proteins_faa):
            if seq_id not in candidate_ids:
                continue
            clean_seq = sanitize_protein_seq(seq)
            out.write(f">{seq_id}\n")
            out.write(_wrap(clean_seq) + "\n")


def load_candidate_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    for idx, raw in enumerate(path.read_text().splitlines()):
        if idx == 0 or not raw:
            continue
        ids.add(raw.split("\t", 1)[0])
    return ids


def iter_fasta(path: Path) -> Iterable[tuple[str, str]]:
    seq_id: Optional[str] = None
    seq: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(seq)
                header = line[1:].strip()
                seq_id = header.split()[0]
                seq = []
                continue
            seq.append(line.strip())
    if seq_id is not None:
        yield seq_id, "".join(seq)


def _wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def sanitize_protein_seq(seq: str) -> str:
    return "".join(ch for ch in seq.upper() if "A" <= ch <= "Z")


def write_structure_map(*, candidates_tsv: Path, structures_dir: Path, out_tsv: Path) -> None:
    if not structures_dir.exists():
        raise MtaseMotifError(f"structures_dir does not exist: {structures_dir}")
    if not structures_dir.is_dir():
        raise MtaseMotifError(f"structures_dir is not a directory: {structures_dir}")

    rows: list[tuple[str, str]] = []
    for candidate_id in sorted(load_candidate_ids(candidates_tsv)):
        structure = find_structure_file(structures_dir, candidate_id)
        rows.append((candidate_id, str(structure) if structure is not None else ""))

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w") as handle:
        handle.write("candidate_id\tstructure_path\n")
        for candidate_id, structure_path in rows:
            handle.write(f"{candidate_id}\t{structure_path}\n")


def find_structure_file(structures_dir: Path, candidate_id: str) -> Optional[Path]:
    for ext in STRUCTURE_EXTS:
        exact = structures_dir / f"{candidate_id}{ext}"
        if exact.exists():
            return exact

    for ext in STRUCTURE_EXTS:
        for path in sorted(structures_dir.glob(f"{candidate_id}*{ext}")):
            if _matches_structure_candidate(candidate_id, path, ext):
                return path
    return None


def run_foldseek_search(
    *,
    structure_map: Path,
    out_tsv: Path,
    foldseek_db: Path,
    threads: int,
) -> None:
    foldseek = shutil.which("foldseek")
    if foldseek is None:
        raise MtaseMotifError("foldseek not found in PATH (required for structure search)")
    header = "query\ttarget\tevalue\tbits\talntmscore\n"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    if not foldseek_db.exists():
        raise MtaseMotifError(f"Foldseek DB not found: {foldseek_db}")

    query_dir = out_tsv.parent / "query_structures"
    if query_dir.exists():
        shutil.rmtree(query_dir)
    query_dir.mkdir(parents=True, exist_ok=True)

    query_count = 0
    for candidate_id, path in load_structure_map(structure_map):
        if not path:
            continue
        src = Path(path)
        if not src.exists():
            continue
        dest = query_dir / _query_filename(candidate_id, src)
        _link_or_copy(src, dest)
        query_count += 1

    if query_count == 0:
        out_tsv.write_text(header)
        return

    tmp_dir = out_tsv.parent / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    raw_out = out_tsv.parent / "foldseek_raw.tsv"

    argv = [
        foldseek,
        "easy-search",
        str(query_dir),
        str(foldseek_db),
        str(raw_out),
        str(tmp_dir),
        "--threads",
        str(threads),
        "--format-output",
        "query,target,evalue,bits,alntmscore",
    ]
    proc = subprocess.run(argv, capture_output=True, text=True)
    if proc.returncode != 0:
        raise MtaseMotifError(proc.stderr.strip() or f"foldseek failed ({proc.returncode})")

    with out_tsv.open("w") as handle:
        handle.write(header)
        if raw_out.exists() and raw_out.stat().st_size > 0:
            handle.write(raw_out.read_text())


def load_structure_map(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    lines = path.read_text().splitlines()
    if not lines:
        return rows
    start = 1 if lines[0].startswith("candidate_id\t") else 0
    for raw in lines[start:]:
        if not raw:
            continue
        parts = raw.split("\t")
        if len(parts) < 2:
            continue
        candidate_id = parts[0].strip()
        structure = parts[1].strip()
        if candidate_id:
            rows.append((candidate_id, structure))
    return rows


def _query_filename(candidate_id: str, src: Path) -> str:
    suffix = "".join(src.suffixes)
    if not suffix:
        suffix = ".pdb"
    return f"{candidate_id}{suffix}"


def _link_or_copy(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        dest.unlink()
    try:
        dest.symlink_to(src)
    except Exception:
        shutil.copy2(src, dest)


def _matches_structure_candidate(candidate_id: str, path: Path, ext: str) -> bool:
    name = path.name
    if not name.endswith(ext):
        return False
    base = name[: -len(ext)]
    if base == candidate_id:
        return True
    if not base.startswith(candidate_id):
        return False
    suffix = base[len(candidate_id) :]
    return bool(suffix) and suffix[0] in {"_", "-", "."}
