from __future__ import annotations

from pathlib import Path

from mtase_motif.mmcif import extract_dna_entity_sequences


def test_extract_dna_entity_sequences_from_entity_poly_loop(tmp_path: Path) -> None:
    cif = tmp_path / "test.cif"
    cif.write_text(
        "data_test\n"
        "#\n"
        "loop_\n"
        "_entity_poly.entity_id\n"
        "_entity_poly.type\n"
        "_entity_poly.pdbx_seq_one_letter_code_can\n"
        "1 'polypeptide(L)'\n"
        ";MSTNPKPQR\n"
        ";\n"
        "2 polydeoxyribonucleotide\n"
        ";ACGTACGT\n"
        ";\n"
        "3 polydeoxyribonucleotide\n"
        ";TTTTCCCC\n"
        ";\n"
        "#\n"
    )
    seqs = extract_dna_entity_sequences(cif)
    assert seqs == ["ACGTACGT", "TTTTCCCC"]


def test_extract_dna_entity_sequences_strips_parenthetical_mods(tmp_path: Path) -> None:
    cif = tmp_path / "mods.cif"
    cif.write_text(
        "data_mods\n"
        "loop_\n"
        "_entity_poly.entity_id\n"
        "_entity_poly.type\n"
        "_entity_poly.pdbx_seq_one_letter_code\n"
        "1 polydeoxyribonucleotide\n"
        ";A C(5MC) G T\n"
        ";\n"
    )
    assert extract_dna_entity_sequences(cif) == ["ACGT"]

