from __future__ import annotations

from mtase_motif.panel import (
    add_weighted_kmers,
    canonical_kmer,
    choose_k_by_enrichment,
    extract_pdb_id,
    reverse_complement,
    select_top_kmers,
)


def test_reverse_complement() -> None:
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAN") == "NTT"


def test_canonical_kmer() -> None:
    assert canonical_kmer("GAATTC") == "GAATTC"
    assert canonical_kmer("GATC") == "GATC"
    assert canonical_kmer("ACG") == "ACG"


def test_add_weighted_kmers_counts_and_canonicalizes() -> None:
    counts: dict[str, float] = {}
    add_weighted_kmers(counts, "ACGTACGT", k=4, weight=2.0)
    assert counts == {"ACGT": 4.0, "CGTA": 4.0, "GTAC": 2.0}


def test_select_top_kmers_stable_order() -> None:
    kmers, weights = select_top_kmers({"AAAA": 2.0, "AAAT": 2.0, "TTTT": 3.0}, limit=2)
    assert kmers == ["TTTT", "AAAA"]
    assert weights == [3.0, 2.0]


def test_choose_k_by_enrichment() -> None:
    counts_by_k = {
        6: {"AAAAAA": 10.0, "CCCCCC": 5.0},
        7: {"AAAAAAA": 5.0, "CCCCCCC": 5.0},
    }
    assert choose_k_by_enrichment(counts_by_k) == 6


def test_extract_pdb_id() -> None:
    assert extract_pdb_id("1abc_A") == "1abc"
    assert extract_pdb_id("/tmp/1ABC.cif.gz") == "1abc"
    assert extract_pdb_id("nope") is None
