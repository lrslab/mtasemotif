from mtase_motif.methylation import methylation_labels, normalize_methylation


def test_normalize_methylation_uses_parenthesized_base_class() -> None:
    assert methylation_labels("5(6)") == {"m6A"}
    assert normalize_methylation("5(6)") == "m6A"
    assert methylation_labels("2(5)") == {"m5C"}
    assert normalize_methylation("2(5)") == "m5C"


def test_normalize_methylation_preserves_canonical_multi_class_labels() -> None:
    assert methylation_labels("m6A/m4C") == {"m6A", "m4C"}
    assert normalize_methylation("m6A/m4C") == "m6A/m4C"
