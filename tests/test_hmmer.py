from __future__ import annotations

from pathlib import Path

from mtase_motif.hmmer import parse_domtblout


def _domtbl_line(model_name: str, model_acc: str, query: str, ievalue: float) -> str:
    cols = [
        model_name,
        model_acc,
        "100",
        query,
        "-",
        "300",
        "1e-30",
        "200.0",
        "0.0",
        "1",
        "1",
        f"{ievalue:g}",
        f"{ievalue:g}",
        "200.0",
        "0.0",
        "1",
        "100",
        "1",
        "100",
        "1",
        "100",
        "0.99",
        "desc",
    ]
    return " ".join(cols) + "\n"


def test_parse_domtblout_basic(tmp_path: Path) -> None:
    domtbl = tmp_path / "hits.domtblout"
    domtbl.write_text(
        "# comment\n"
        + _domtbl_line("N6_N4_Mtase", "PF01555.1", "gene1", 1e-50)
        + _domtbl_line("SomeOther", "PF00000.1", "gene1", 1e-2)
    )

    hits = parse_domtblout(domtbl, 1e-5, source="pfam")
    assert "gene1" in hits
    assert len(hits["gene1"]) == 1
    hit = hits["gene1"][0]
    assert hit.model_acc == "PF01555"
    assert hit.model_name == "N6_N4_Mtase"
    assert hit.source == "pfam"


def test_parse_domtblout_hmmsearch_uses_target_as_query_id(tmp_path: Path) -> None:
    domtbl = tmp_path / "hmmsearch.domtblout"
    domtbl.write_text(
        "protein123 - 300 TIGR01234 TIGR01234 100 1e-30 200.0 0.0 1 1 1e-40 1e-40 200.0 0.0 1 100 1 100 1 100 0.99 desc\n"
    )

    hits = parse_domtblout(domtbl, 1e-5, source="tigrfams", program="hmmsearch")
    assert "protein123" in hits
    assert len(hits["protein123"]) == 1
    hit = hits["protein123"][0]
    assert hit.model_acc == "TIGR01234"
    assert hit.model_name == "TIGR01234"
    assert hit.source == "tigrfams"
