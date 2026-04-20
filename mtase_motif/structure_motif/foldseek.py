from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class FoldseekHit:
    target: str
    evalue: float
    bits: float
    alntmscore: float


def confidence_from_tmscore(tmscore: float) -> str:
    if tmscore >= 0.7:
        return "high"
    if tmscore >= 0.55:
        return "medium"
    return "low"


def choose_best_foldseek_hit(hits: list[FoldseekHit], *, min_alntmscore: float) -> FoldseekHit | None:
    passing = [h for h in hits if h.alntmscore >= min_alntmscore]
    if not passing:
        return None
    passing.sort(key=lambda h: (-h.alntmscore, h.evalue, -h.bits, h.target))
    return passing[0]


def load_foldseek_hits(path: Path) -> dict[str, list[FoldseekHit]]:
    hits: dict[str, list[FoldseekHit]] = {}
    lines = path.read_text().splitlines()
    if not lines:
        return hits

    start = 1 if lines[0].startswith("query\t") else 0
    for raw in lines[start:]:
        if not raw.strip():
            continue
        parts = raw.split("\t")
        if len(parts) < 5:
            continue
        query_raw, target = parts[0], parts[1]
        try:
            evalue = float(parts[2])
            bits = float(parts[3])
            alntmscore = float(parts[4])
        except ValueError:
            continue
        query = _query_to_candidate_id(query_raw)
        if not query:
            continue
        hits.setdefault(query, []).append(
            FoldseekHit(target=target, evalue=evalue, bits=bits, alntmscore=alntmscore)
        )
    return hits


def load_foldseek_labels(path: Path) -> dict[str, tuple[str, str]]:
    labels: dict[str, tuple[str, str]] = {}
    for i, raw in enumerate(path.read_text().splitlines()):
        if i == 0 and raw.startswith("target"):
            continue
        if not raw.strip() or raw.startswith("#"):
            continue
        parts = raw.split("\t")
        if len(parts) < 2:
            continue
        target = parts[0].strip()
        motif_iupac = parts[1].strip()
        methylation = parts[2].strip() if len(parts) > 2 else ""
        if target and motif_iupac:
            labels[target] = (motif_iupac, methylation)
    return labels


def _query_to_candidate_id(query: str) -> str:
    name = Path(query).name
    for ext in (".gz", ".pdb", ".cif", ".mmcif"):
        if name.endswith(ext):
            name = name[: -len(ext)]
    for ext in (".pdb", ".cif", ".mmcif"):
        if name.endswith(ext):
            name = name[: -len(ext)]
    return name
