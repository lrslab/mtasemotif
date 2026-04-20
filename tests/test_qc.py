from __future__ import annotations

from pathlib import Path

from mtase_motif.qc import can_use_exact_iupac_fallback, parse_fimo_hits, write_exact_iupac_fimo


def test_parse_fimo_hits_counts_and_balance(tmp_path: Path) -> None:
    fimo = tmp_path / "fimo.tsv"
    fimo.write_text(
        "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n"
        "M1\t.\tctg1\t1\t4\t+\t10\t1e-6\t1e-6\tGATC\n"
        "M1\t.\tctg1\t10\t13\t-\t10\t1e-6\t1e-6\tGATC\n"
    )
    total, density, balance = parse_fimo_hits(fimo, genome_len=1_000_000)
    assert total == 2
    assert abs(density - 2.0) < 1e-9
    assert balance == 0.0


def test_parse_fimo_hits_skips_comments_and_deduplicates_palindromic_sites(tmp_path: Path) -> None:
    fimo = tmp_path / "fimo.tsv"
    fimo.write_text(
        "# comment\n"
        "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n"
        "M1\t.\tctg1\t1\t5\t+\t10\t1e-6\t1e-6\tCCAGG\n"
        "M1\t.\tctg1\t1\t5\t-\t10\t1e-6\t1e-6\tCCTGG\n"
    )

    total, density, balance = parse_fimo_hits(fimo, genome_len=1_000_000)
    assert total == 1
    assert abs(density - 1.0) < 1e-9
    assert balance == 0.0


def test_write_exact_iupac_fimo_emits_expected_hits(tmp_path: Path) -> None:
    genome = tmp_path / "genome.fa"
    fimo = tmp_path / "fimo.tsv"
    genome.write_text(">ctg1\nACCAGGCCTGG\n")

    write_exact_iupac_fimo(fimo, genome_fa=genome, motif_id="cand1", motif_iupac="CCWGG")

    lines = [line for line in fimo.read_text().splitlines() if line and not line.startswith("#")]
    assert lines[0].startswith("motif_id\tmotif_alt_id\tsequence_name")
    assert len(lines) == 5
    total, density, balance = parse_fimo_hits(fimo, genome_len=1_000_000)
    assert total == 2
    assert abs(density - 2.0) < 1e-9
    assert balance == 0.0


def test_can_use_exact_iupac_fallback_guards_highly_degenerate_motifs() -> None:
    assert can_use_exact_iupac_fallback("CCWGG")
    assert can_use_exact_iupac_fallback("NNNNNNNN") is False
