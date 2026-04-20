from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from mtase_motif.config import DbLayout
from mtase_motif.db.providers.rebase import RebaseProvider
import mtase_motif.db.providers.rebase as rebase_provider
from mtase_motif.util import MtaseMotifError


def _write_rebase_emboss_files(src: Path) -> None:
    (src / "emboss_e.txt").write_text("EcoRI GAATTC\n")
    (src / "emboss_r.txt").write_text(
        "\n".join(
            [
                "EcoRI",
                "Escherichia coli",
                "",
                "5(6)",
                "source",
                "suppliers",
                "0",
                "//",
                "",
            ]
        )
    )
    (src / "emboss_s.txt").write_text("# placeholder\n")


def test_rebase_provider_fetch_and_index_smoke(tmp_path: Path) -> None:
    if shutil.which("mmseqs") is None and shutil.which("makeblastdb") is None:
        pytest.skip("mmseqs or makeblastdb required for indexing")

    src = tmp_path / "rebase_src"
    src.mkdir()
    _write_rebase_emboss_files(src)
    (src / "rebase_proteins.faa").write_text(
        ">EcoRI\nMSTNPKPQRK\n>EcoRI_dup\nMSTNPKPQRK\n>Other\nMAAA\n"
    )

    layout = DbLayout(db_dir=tmp_path / "db")
    layout.ensure_dirs()
    provider = RebaseProvider(layout)
    pm = provider.fetch(source=src)
    provider.index()

    versions = sorted([p for p in layout.rebase_dir.iterdir() if p.is_dir()])
    assert versions, "expected at least one rebase/<version> directory"
    latest = versions[-1]
    assert (latest / "rebase_enzymes.tsv").exists()
    assert (latest / "rebase_proteins.faa").exists()
    assert (latest / "rebase_protein_map.tsv").exists()
    assert (latest / "rebase_motif_proteins.faa").exists()
    assert (latest / "rebase_cluster_labels.tsv").exists()
    assert (latest / "manifest.json").exists()
    protein_map_lines = (latest / "rebase_protein_map.tsv").read_text().splitlines()
    assert protein_map_lines[0].startswith("protein_id\tenzyme_id\tmatched_name")
    assert any(line.startswith("EcoRI\tEcoRI\t") for line in protein_map_lines[1:])
    enzyme_lines = (latest / "rebase_enzymes.tsv").read_text().splitlines()
    enzyme_header = enzyme_lines[0].split("\t")
    enzyme_row = dict(zip(enzyme_header, enzyme_lines[1].split("\t")))
    assert enzyme_row["methylation"] == "m6A"
    assert enzyme_row["methylation_raw"] == "5(6)"
    motif_proteins = (latest / "rebase_motif_proteins.faa").read_text()
    assert ">EcoRI\n" in motif_proteins
    cluster_lines = (latest / "rebase_cluster_labels.tsv").read_text().splitlines()
    assert cluster_lines[0].startswith("protein_id\tcluster_id\trepresentative_id")
    assert any("\tGAATTC\tm6A\thigh\tdirect_duplicate_consensus\t" in line for line in cluster_lines[1:])
    assert any(url.endswith("/rebase_proteins.faa") for url in pm.urls)
    assert "Normalized proteins from rebase_proteins.faa" in (pm.notes or "")


def test_download_remote_proteins_merges_multiple_sets(monkeypatch, tmp_path: Path) -> None:
    urls = (
        "ftp://ftp.neb.com/pub/rebase/All_REBASE_Gold_Standards_Protein.txt",
        "ftp://ftp.neb.com/pub/rebase/Other_methyltransferase_genes_Protein.txt",
    )
    payloads = {
        urls[0]: ">M.A\nMPEPTIDE\n>M.B\nMSEQ\n",
        urls[1]: ">M.B duplicate\nMSEQ\n>M.C\nMNEWSEQ\n",
    }

    def fake_download_url(url: str, dest: Path) -> None:
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(payloads[url])

    monkeypatch.setattr(rebase_provider, "REMOTE_PROTEIN_URLS", urls)
    monkeypatch.setattr(rebase_provider, "_download_url", fake_download_url)

    used_urls, note = rebase_provider._download_remote_proteins(tmp_path)

    proteins = tmp_path / "rebase_proteins.faa"
    text = proteins.read_text()
    assert used_urls == list(urls)
    assert note is not None and "3 unique sequences" in note
    assert text.count("\n>") == 2
    assert ">M.A\n" in text
    assert ">M.B\n" in text
    assert ">M.C\n" in text


def test_rebase_provider_fetch_merges_multiple_local_protein_files(tmp_path: Path) -> None:
    src = tmp_path / "rebase_src"
    src.mkdir()
    _write_rebase_emboss_files(src)
    (src / "Type_II_methyltransferase_genes_Protein.txt").write_text(
        ">M.One\nMPEPTIDE\n>M.Two\nMSECOND\n"
    )
    (src / "Other_methyltransferase_genes_Protein.txt").write_text(
        ">M.One duplicate\nMPEPTIDE\n>M.Three\nMTHIRD\n"
    )

    layout = DbLayout(db_dir=tmp_path / "db")
    layout.ensure_dirs()
    provider = RebaseProvider(layout)
    pm = provider.fetch(source=src)

    versions = sorted([p for p in layout.rebase_dir.iterdir() if p.is_dir()])
    latest = versions[-1]
    proteins = (latest / "rebase_proteins.faa").read_text().splitlines()
    headers = [line for line in proteins if line.startswith(">")]

    assert len(headers) == 3
    assert headers == [">M.One", ">M.Two", ">M.Three"]
    assert any(url.endswith("Type_II_methyltransferase_genes_Protein.txt") for url in pm.urls)
    assert any(url.endswith("Other_methyltransferase_genes_Protein.txt") for url in pm.urls)
    assert "Merged 2 REBASE protein source files into 3 unique sequences" in (pm.notes or "")


def test_rebase_provider_prefers_canonical_protein_fasta_when_present(tmp_path: Path) -> None:
    src = tmp_path / "rebase_src"
    src.mkdir()
    _write_rebase_emboss_files(src)
    (src / "rebase_proteins.faa").write_text(">Canonical\nMPEPTIDE\n")
    (src / "Type_II_methyltransferase_genes_Protein.txt").write_text(">Ignored\nMOTHER\n")

    layout = DbLayout(db_dir=tmp_path / "db")
    layout.ensure_dirs()
    provider = RebaseProvider(layout)
    pm = provider.fetch(source=src)

    versions = sorted([p for p in layout.rebase_dir.iterdir() if p.is_dir()])
    latest = versions[-1]
    text = (latest / "rebase_proteins.faa").read_text()

    assert text == ">Canonical\nMPEPTIDE\n"
    assert any(url.endswith("/rebase_proteins.faa") for url in pm.urls)
    assert not any(url.endswith("Type_II_methyltransferase_genes_Protein.txt") for url in pm.urls)
    assert "Normalized proteins from rebase_proteins.faa" in (pm.notes or "")


def test_pick_emboss_files_rejects_mixed_release_suffixes(tmp_path: Path) -> None:
    src = tmp_path / "rebase_src"
    src.mkdir()
    (src / "emboss_e.001").write_text("EcoRI GAATTC\n")
    (src / "emboss_r.002").write_text("EcoRI\n")
    (src / "emboss_s.001").write_text("# placeholder\n")

    with pytest.raises(MtaseMotifError, match="coherent release"):
        rebase_provider._pick_emboss_files(src)


def test_pick_emboss_files_prefers_latest_complete_numeric_release(tmp_path: Path) -> None:
    src = tmp_path / "rebase_src"
    src.mkdir()
    for suffix in ("001", "002"):
        (src / f"emboss_e.{suffix}").write_text("EcoRI GAATTC\n")
        (src / f"emboss_r.{suffix}").write_text("EcoRI\n")
        (src / f"emboss_s.{suffix}").write_text("# placeholder\n")

    emboss_e, emboss_r, emboss_s, version = rebase_provider._pick_emboss_files(src)

    assert emboss_e.name == "emboss_e.002"
    assert emboss_r.name == "emboss_r.002"
    assert emboss_s.name == "emboss_s.002"
    assert version == "002"
