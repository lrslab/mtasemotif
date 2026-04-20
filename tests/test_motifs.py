from __future__ import annotations

from mtase_motif.motifs import (
    canonicalize_iupac_orientation,
    classify_iupac_motif,
    infer_mod_position,
    iupac_degeneracy,
    is_palindromic_iupac,
    iupac_to_pwm,
    matches_iupac_sequence,
    normalize_iupac,
    pwm_to_iupac,
    reverse_complement_iupac,
    sequences_to_pwm,
)


def test_normalize_iupac_strips_emboss_markers() -> None:
    assert normalize_iupac("G^AATTC") == "GAATTC"
    assert normalize_iupac("CCWGG") == "CCWGG"


def test_iupac_to_pwm_rows_sum_to_one() -> None:
    pwm = iupac_to_pwm("AN")
    assert len(pwm) == 2
    for row in pwm:
        assert abs(sum(row) - 1.0) < 1e-9


def test_is_palindromic_iupac() -> None:
    assert is_palindromic_iupac("GAATTC")
    assert is_palindromic_iupac("GATC")
    assert is_palindromic_iupac("CCWGG")
    assert not is_palindromic_iupac("GACT")


def test_sequences_to_pwm_and_pwm_to_iupac_roundtrip() -> None:
    pwm = sequences_to_pwm(["AAAA", "AAAT", "AAAT"])
    assert len(pwm) == 4
    for row in pwm:
        assert abs(sum(row) - 1.0) < 1e-9
    assert pwm_to_iupac(pwm) == "AAAW"


def test_reverse_complement_iupac_and_matching() -> None:
    assert reverse_complement_iupac("CCWGG") == "CCWGG"
    assert reverse_complement_iupac("GANTC") == "GANTC"
    assert reverse_complement_iupac("CATG") == "CATG"
    assert matches_iupac_sequence("CCAGG", "CCWGG")
    assert matches_iupac_sequence("CCTGG", "CCWGG")
    assert not matches_iupac_sequence("CCCGG", "CCWGG")


def test_iupac_degeneracy_counts_expansions() -> None:
    assert iupac_degeneracy("GATC") == 1
    assert iupac_degeneracy("CCWGG") == 2
    assert iupac_degeneracy("NNNN") == 256


def test_canonicalize_iupac_orientation_picks_stable_strand() -> None:
    assert canonicalize_iupac_orientation("CCWGG") == "CCWGG"
    assert canonicalize_iupac_orientation("TGTA") == "TACA"


def test_classify_iupac_motif_distinguishes_bipartite_and_palindromic() -> None:
    assert classify_iupac_motif("CCWGG") == "palindromic"
    assert classify_iupac_motif("GACT") == "non_palindromic"
    assert classify_iupac_motif("AACNNNNGTGC") == "bipartite"


def test_infer_mod_position_only_returns_unique_supported_site() -> None:
    assert infer_mod_position("GATC", "m6A") == 2
    assert infer_mod_position("CCWGG", "m5C") is None
    assert infer_mod_position("CCSGG", "m4C") is None
    assert infer_mod_position("GATC", "m6A/m4C") is None
