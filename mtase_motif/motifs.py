from __future__ import annotations

import re
from pathlib import Path

from mtase_motif.methylation import methylation_labels


IUPAC_TO_BASES: dict[str, set[str]] = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

_COMPLEMENT: dict[str, str] = {"A": "T", "C": "G", "G": "C", "T": "A"}
_BASES: tuple[str, ...] = ("A", "C", "G", "T")

_BASESET_TO_IUPAC: dict[frozenset[str], str] = {}
for _sym in ("A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"):
    _BASESET_TO_IUPAC[frozenset(IUPAC_TO_BASES[_sym])] = _sym


def normalize_iupac(seq: str) -> str:
    out: list[str] = []
    for ch in seq.upper():
        if ch in IUPAC_TO_BASES:
            out.append(ch)
        elif ch in {"^", "/", "(", ")", "[", "]", "{", "}", "-", "_"}:
            continue
        elif ch.isspace():
            continue
        else:
            out.append("N")
    return "".join(out)


def iupac_to_pwm(seq: str) -> list[list[float]]:
    matrix: list[list[float]] = []
    for ch in seq:
        allowed = IUPAC_TO_BASES.get(ch, IUPAC_TO_BASES["N"])
        p = 1.0 / float(len(allowed))
        matrix.append([p if b in allowed else 0.0 for b in _BASES])
    return matrix


def iupac_degeneracy(seq: str) -> int:
    total = 1
    for ch in normalize_iupac(seq):
        total *= len(IUPAC_TO_BASES.get(ch, IUPAC_TO_BASES["N"]))
    return total


def reverse_complement_iupac(seq: str) -> str:
    out: list[str] = []
    for ch in reversed(normalize_iupac(seq)):
        bases = IUPAC_TO_BASES.get(ch, IUPAC_TO_BASES["N"])
        out.append(_BASESET_TO_IUPAC.get(frozenset(_complement_bases(bases)), "N"))
    return "".join(out)


def canonicalize_iupac_orientation(seq: str) -> str:
    norm = normalize_iupac(seq)
    if not norm:
        return ""
    rc = reverse_complement_iupac(norm)
    return min(norm, rc)


def matches_iupac_sequence(seq: str, motif_iupac: str) -> bool:
    seq_norm = seq.upper().replace("U", "T")
    motif_norm = normalize_iupac(motif_iupac)
    if len(seq_norm) != len(motif_norm):
        return False
    for base, sym in zip(seq_norm, motif_norm):
        if base not in IUPAC_TO_BASES.get(sym, IUPAC_TO_BASES["N"]):
            return False
    return True


def write_meme_pwm(path: Path, *, motif_id: str, motif_iupac: str) -> None:
    matrix = iupac_to_pwm(motif_iupac)
    write_meme_matrix(path, motif_id=motif_id, matrix=matrix, nsites=20)


def write_meme_matrix(path: Path, *, motif_id: str, matrix: list[list[float]], nsites: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("MEME version 4\n\n")
        f.write("ALPHABET= ACGT\n\n")
        f.write("strands: + -\n\n")
        f.write("Background letter frequencies (from uniform background):\n")
        f.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")
        f.write(f"MOTIF {motif_id}\n")
        f.write(
            f"letter-probability matrix: alength= 4 w= {len(matrix)} nsites= {nsites} E= 0\n"
        )
        for row in matrix:
            f.write(" ".join(f"{v:.6f}" for v in row) + "\n")


def sequences_to_pwm(
    sequences: list[str],
    *,
    weights: list[float] | None = None,
    pseudocount: float = 0.01,
) -> list[list[float]]:
    if not sequences:
        raise ValueError("sequences_to_pwm requires at least one sequence")
    width = len(sequences[0])
    if width == 0:
        raise ValueError("sequences_to_pwm requires non-empty sequences")
    if any(len(s) != width for s in sequences):
        raise ValueError("sequences_to_pwm requires all sequences to have the same length")
    if weights is not None and len(weights) != len(sequences):
        raise ValueError("weights length must match sequences length")

    matrix: list[list[float]] = []
    for pos in range(width):
        counts = {b: float(pseudocount) for b in _BASES}
        for i, seq in enumerate(sequences):
            base = seq[pos].upper()
            if base == "U":
                base = "T"
            if base not in counts:
                continue
            w = float(weights[i]) if weights is not None else 1.0
            counts[base] += w
        total = sum(counts.values())
        matrix.append([counts[b] / total for b in _BASES])
    return matrix


def pwm_to_iupac(matrix: list[list[float]], *, min_prob: float = 0.2) -> str:
    if not matrix:
        return ""
    out: list[str] = []
    for row in matrix:
        if len(row) != 4:
            raise ValueError("pwm_to_iupac expects rows of length 4 (A,C,G,T)")
        bases = {b for b, p in zip(_BASES, row) if p >= min_prob}
        if not bases:
            bases = {_BASES[max(range(4), key=lambda i: row[i])]}
        out.append(_BASESET_TO_IUPAC.get(frozenset(bases), "N"))
    return "".join(out)


def is_palindromic_iupac(seq: str) -> bool:
    if not seq:
        return False
    if len(seq) == 1:
        return True

    bases = [IUPAC_TO_BASES.get(ch, IUPAC_TO_BASES["N"]) for ch in seq.upper()]
    for i in range(len(bases)):
        left = bases[i]
        right = bases[len(bases) - 1 - i]
        if left != _complement_bases(right):
            return False
    return True


def classify_iupac_motif(seq: str) -> str:
    norm = normalize_iupac(seq)
    if not norm:
        return ""
    if re.search(r"[^N]N{2,}[^N]", norm):
        return "bipartite"
    if is_palindromic_iupac(norm):
        return "palindromic"
    return "non_palindromic"


def infer_mod_position(seq: str, methylation: str) -> int | None:
    norm = normalize_iupac(seq)
    labels = methylation_labels(methylation)
    if not norm or len(labels) != 1:
        return None

    label = next(iter(labels))
    target_base = "A" if label == "m6A" else "C"
    matching_positions = [
        idx + 1 for idx, sym in enumerate(norm) if target_base in IUPAC_TO_BASES.get(sym, IUPAC_TO_BASES["N"])
    ]
    if len(matching_positions) != 1:
        return None
    return matching_positions[0]


def _complement_bases(bases: set[str]) -> set[str]:
    return {_COMPLEMENT.get(b, "N") for b in bases}
