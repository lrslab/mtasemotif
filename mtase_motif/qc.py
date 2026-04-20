from __future__ import annotations

from pathlib import Path

from mtase_motif.motifs import (
    iupac_degeneracy,
    matches_iupac_sequence,
    normalize_iupac,
    reverse_complement_iupac,
)
from mtase_motif.panel import reverse_complement


def can_use_exact_iupac_fallback(motif_iupac: str, *, max_degeneracy: int = 256) -> bool:
    motif = normalize_iupac(motif_iupac)
    return bool(motif) and iupac_degeneracy(motif) <= max_degeneracy


def write_exact_iupac_fimo(
    fimo_tsv: Path,
    *,
    genome_fa: Path,
    motif_id: str,
    motif_iupac: str,
) -> None:
    motif = normalize_iupac(motif_iupac)
    rc_motif = reverse_complement_iupac(motif)

    fimo_tsv.parent.mkdir(parents=True, exist_ok=True)
    with fimo_tsv.open("w") as out:
        out.write("# mtase-motif exact IUPAC fallback because FIMO returned no hits\n")
        out.write(
            "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n"
        )
        for seq_name, seq in _iter_fasta(genome_fa):
            width = len(motif)
            if width == 0 or len(seq) < width:
                continue
            for i in range(len(seq) - width + 1):
                window = seq[i : i + width]
                plus = matches_iupac_sequence(window, motif)
                minus = matches_iupac_sequence(window, rc_motif)
                start = str(i + 1)
                stop = str(i + width)
                if plus:
                    out.write(
                        f"{motif_id}\t\t{seq_name}\t{start}\t{stop}\t+\t0\t1\t1\t{window}\n"
                    )
                if minus:
                    out.write(
                        f"{motif_id}\t\t{seq_name}\t{start}\t{stop}\t-\t0\t1\t1\t{reverse_complement(window)}\n"
                    )


def _iter_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = ""
    seq_chunks: list[str] = []
    with path.open() as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    records.append((name, "".join(seq_chunks).upper().replace("U", "T")))
                name = line[1:].split()[0]
                seq_chunks = []
                continue
            seq_chunks.append(line)
    if name:
        records.append((name, "".join(seq_chunks).upper().replace("U", "T")))
    return records


def parse_fimo_hits(fimo_tsv: Path, genome_len: int) -> tuple[int, float, float]:
    if not fimo_tsv.exists():
        return 0, 0.0, 0.0
    plus = 0
    minus = 0
    unique_sites: set[tuple[str, str, str, str]] = set()
    for raw in fimo_tsv.read_text().splitlines():
        if not raw or raw.startswith("#"):
            continue
        parts = raw.split("\t")
        if parts[0] == "motif_id":
            continue
        if len(parts) < 7:
            continue
        unique_sites.add((parts[0], parts[2], parts[3], parts[4]))
        strand = parts[5]
        if strand == "+":
            plus += 1
        elif strand == "-":
            minus += 1

    total = len(unique_sites)
    density = 0.0 if genome_len <= 0 else total / (genome_len / 1e6)
    strand_total = plus + minus
    strand_balance = 0.0 if strand_total == 0 else abs(plus - minus) / float(strand_total)
    return total, density, strand_balance
