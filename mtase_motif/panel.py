from __future__ import annotations

import re
from typing import Dict, Optional, Tuple


_RC: dict[str, str] = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

_PDB_ID_RE = re.compile(r"(?i)(?<![A-Z0-9])[0-9][A-Z0-9]{3}(?![A-Z0-9])")


def reverse_complement(seq: str) -> str:
    s = seq.upper()
    return "".join(_RC.get(ch, "N") for ch in reversed(s))


def canonical_kmer(kmer: str) -> str:
    k = kmer.upper()
    rc = reverse_complement(k)
    return rc if rc < k else k


def add_weighted_kmers(
    counts: Dict[str, float],
    seq: str,
    *,
    k: int,
    weight: float,
    canonicalize: bool = True,
) -> None:
    if k <= 0:
        return
    s = seq.upper()
    if len(s) < k:
        return
    for i in range(len(s) - k + 1):
        kmer = s[i : i + k]
        if "N" in kmer:
            continue
        if canonicalize:
            kmer = canonical_kmer(kmer)
        counts[kmer] = counts.get(kmer, 0.0) + float(weight)


def select_top_kmers(counts: Dict[str, float], *, limit: int) -> Tuple[list[str], list[float]]:
    if limit <= 0 or not counts:
        return [], []
    items = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    kmers = [k for k, _w in items[:limit]]
    weights = [float(w) for _k, w in items[:limit]]
    return kmers, weights


def choose_k_by_enrichment(counts_by_k: Dict[int, Dict[str, float]]) -> Optional[int]:
    best_k: Optional[int] = None
    best_score: tuple[float, float, int] | None = None
    for k, counts in counts_by_k.items():
        if not counts:
            continue
        total = float(sum(counts.values()))
        if total <= 0:
            continue
        top = float(max(counts.values()))
        score = (top / total, top, -int(k))
        if best_score is None or score > best_score:
            best_score = score
            best_k = int(k)
    return best_k


def extract_pdb_id(value: str) -> Optional[str]:
    if not value:
        return None
    m = _PDB_ID_RE.search(value)
    if not m:
        return None
    return m.group(0).lower()
