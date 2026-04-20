from __future__ import annotations

import json
import shutil
from pathlib import Path
import subprocess
from typing import Optional

import pytest

from mtase_motif.candidate_discovery import build_candidates
from mtase_motif.motif_inference import infer_motifs
from mtase_motif.pipeline_support import (
    extract_candidate_proteins,
    find_structure_file,
    run_foldseek_search,
)
from mtase_motif.rebase_search import run_blast, run_mmseqs, search_rebase
from mtase_motif.scanning import scan_and_summarize
from mtase_motif.structure_motif.template_panel import _download_url
from mtase_motif.util import MtaseMotifError


MMCIF_WITH_DNA = """data_test
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_seq_one_letter_code_can
1 polydeoxyribonucleotide GATC
#
"""


def _candidate_tsv(path: Path) -> None:
    path.write_text(
        "candidate_id\tcontig\tstart\tend\tstrand\tdomains\ttype_hint\tmethylation_hint\tneighborhood_features\tconfidence\trationale\n"
        "cand1\tctg1\t1\t100\t+\tpfam:PF01555(N6_N4_Mtase):1e-40\ttype_II_or_orphan_like\tm6A/m4C\twindow_bp=10000\t0.9\tpfam_hit\n"
    )


def _enzymes_tsv(path: Path) -> None:
    path.write_text(
        "enzyme_id\tnames\ttype\trecognition_seq_iupac\tmethylation\tevidence\n"
        "EcoRI\tEcoRI\tunknown\tGAATTC\tm6A\trebase_emboss\n"
    )


def test_build_candidates_maps_prodigal_protein_ids_back_to_gff_coords(tmp_path: Path) -> None:
    pfam_domtbl = tmp_path / "pfam.domtblout"
    tigr_domtbl = tmp_path / "tigr.domtblout"
    gff = tmp_path / "genes.gff"
    out_tsv = tmp_path / "mtase_candidates.tsv"

    pfam_domtbl.write_text(
        "N6_N4_Mtase PF01555.25 100 ctg1_42 - 300 1e-30 200.0 0.0 1 1 1e-50 1e-50 200.0 0.0 1 100 1 100 1 100 0.99 DNA methylase\n"
    )
    tigr_domtbl.write_text("")
    gff.write_text(
        "ctg1\tProdigal:2.6\tCDS\t100\t500\t.\t+\t0\tID=1_42;partial=00;start_type=ATG;\n"
    )

    build_candidates(
        pfam_domtbl_path=pfam_domtbl,
        tigr_domtbl_path=tigr_domtbl,
        gff_path=gff,
        out_path=out_tsv,
        pfam_evalue_thresh=1e-5,
        tigr_evalue_thresh=1e-5,
        max_domains=5,
    )

    lines = out_tsv.read_text().splitlines()
    assert len(lines) == 2
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["candidate_id"] == "ctg1_42"
    assert row["contig"] == "ctg1"
    assert row["type_hint"] == "type_II_or_orphan_like"


def test_build_candidates_adds_named_rebase_rescue_candidate(
    tmp_path: Path,
    monkeypatch,
) -> None:
    pfam_domtbl = tmp_path / "pfam.domtblout"
    tigr_domtbl = tmp_path / "tigr.domtblout"
    gff = tmp_path / "genes.gff"
    proteins = tmp_path / "proteins.faa"
    rebase_proteins = tmp_path / "rebase_proteins.faa"
    out_tsv = tmp_path / "mtase_candidates.tsv"

    pfam_domtbl.write_text("")
    tigr_domtbl.write_text("")
    gff.write_text(
        "ctg1\tProdigal:2.6\tCDS\t100\t936\t.\t+\t0\tID=1_3567;partial=00;start_type=ATG;\n"
    )
    proteins.write_text(">ctg1_3567\nMPEPTIDE\n")
    rebase_proteins.write_text(">M.EcoKDam\nMPEPTIDE\n>M.EcoPE144Dcm3P\nMOTHER\n")

    def fake_run(argv, capture_output, text):
        assert "-subject" in argv
        return subprocess.CompletedProcess(
            argv,
            0,
            stdout="ctg1_3567\tM.EcoKDam\t0.0\t579\t99.640\t100\t279\t278\n",
            stderr="",
        )

    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/blastp" if name == "blastp" else None)
    monkeypatch.setattr(subprocess, "run", fake_run)

    build_candidates(
        pfam_domtbl_path=pfam_domtbl,
        tigr_domtbl_path=tigr_domtbl,
        gff_path=gff,
        out_path=out_tsv,
        proteins_faa=proteins,
        rebase_proteins=rebase_proteins,
        pfam_evalue_thresh=1e-5,
        tigr_evalue_thresh=1e-5,
        max_domains=5,
        rescue_max_evalue=1e-20,
        rescue_min_pident=70.0,
        rescue_min_qcov=80.0,
        rescue_min_tcov=80.0,
    )

    lines = out_tsv.read_text().splitlines()
    assert len(lines) == 2
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["candidate_id"] == "ctg1_3567"
    assert row["methylation_hint"] == "m6A"
    assert row["domains"] == "rebase_named:M.EcoKDam:0"
    assert "rebase_named_rescue M.EcoKDam" in row["rationale"]


def test_build_candidates_keeps_protein_hits_without_gff_coords(tmp_path: Path) -> None:
    pfam_domtbl = tmp_path / "pfam.domtblout"
    tigr_domtbl = tmp_path / "tigr.domtblout"
    gff = tmp_path / "genes.gff"
    out_tsv = tmp_path / "mtase_candidates.tsv"

    pfam_domtbl.write_text(
        "N6_N4_Mtase PF01555.25 100 prot1 - 300 1e-30 200.0 0.0 1 1 1e-50 1e-50 200.0 0.0 1 100 1 100 1 100 0.99 DNA methylase\n"
    )
    tigr_domtbl.write_text("")
    gff.write_text("# No gene calls produced (proteins were provided via --proteins)\n")

    build_candidates(
        pfam_domtbl_path=pfam_domtbl,
        tigr_domtbl_path=tigr_domtbl,
        gff_path=gff,
        out_path=out_tsv,
        pfam_evalue_thresh=1e-5,
        tigr_evalue_thresh=1e-5,
        max_domains=5,
    )

    lines = out_tsv.read_text().splitlines()
    assert len(lines) == 2
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["candidate_id"] == "prot1"
    assert row["contig"] == ""
    assert row["start"] == ""
    assert row["end"] == ""
    assert row["strand"] == ""
    assert row["type_hint"] == "type_II_or_orphan_like"
    assert "coord_source=unavailable" in row["neighborhood_features"]


def test_build_candidates_uses_hmmsearch_orientation_for_tigrfams(tmp_path: Path) -> None:
    pfam_domtbl = tmp_path / "pfam.domtblout"
    tigr_domtbl = tmp_path / "tigr.domtblout"
    gff = tmp_path / "genes.gff"
    out_tsv = tmp_path / "mtase_candidates.tsv"

    pfam_domtbl.write_text("")
    tigr_domtbl.write_text(
        "ctg1_9 - 300 Dam_methyltransferase TIGR04167 100 1e-30 200.0 0.0 1 1 1e-50 1e-50 200.0 0.0 1 100 1 100 1 100 0.99 Dam methyltransferase\n"
    )
    gff.write_text(
        "ctg1\tProdigal:2.6\tCDS\t100\t900\t.\t+\t0\tID=1_9;partial=00;start_type=ATG;\n"
    )

    build_candidates(
        pfam_domtbl_path=pfam_domtbl,
        tigr_domtbl_path=tigr_domtbl,
        gff_path=gff,
        out_path=out_tsv,
        pfam_evalue_thresh=1e-5,
        tigr_evalue_thresh=1e-5,
        max_domains=5,
    )

    lines = out_tsv.read_text().splitlines()
    assert len(lines) == 2
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["candidate_id"] == "ctg1_9"
    assert row["contig"] == "ctg1"
    assert row["methylation_hint"] == "m6A"


def test_extract_candidate_proteins_strips_terminal_stop_codons(tmp_path: Path) -> None:
    proteins = tmp_path / "proteins.faa"
    candidates = tmp_path / "candidates.tsv"
    out_faa = tmp_path / "mtase_candidates.faa"

    proteins.write_text(">cand1\nMPEPTIDE*\n>other\nMOTHER*\n")
    candidates.write_text(
        "candidate_id\tcontig\tstart\tend\tstrand\tdomains\ttype_hint\tmethylation_hint\tneighborhood_features\tconfidence\trationale\n"
        "cand1\tctg1\t1\t10\t+\tpfam:PF01555(N6_N4_Mtase):1e-50\ttype_II_or_orphan_like\tm6A/m4C\twindow_bp=10000\t1.0\tpfam_hit\n"
    )

    extract_candidate_proteins(proteins_faa=proteins, candidates_tsv=candidates, out_faa=out_faa)

    text = out_faa.read_text()
    assert ">cand1\nMPEPTIDE\n" == text


def test_rebase_search_parses_blast_output_with_length_based_target_coverage(
    tmp_path: Path,
    monkeypatch,
) -> None:
    query = tmp_path / "candidates.faa"
    out_tsv = tmp_path / "rebase_search.tsv"
    blast_prefix = tmp_path / "blast" / "rebase_proteins"

    query.write_text(">cand1\nMPEPTIDE\n")
    blast_prefix.parent.mkdir(parents=True)
    for ext in (".pin", ".psq", ".phr"):
        Path(str(blast_prefix) + ext).write_text("")

    def fake_run(argv, capture_output, text):
        assert "qcovs slen length" in argv[-1]
        return subprocess.CompletedProcess(
            argv,
            0,
            stdout="cand1\tM.EcoDcm\t0.0\t987\t100.000\t100\t472\t472\n",
            stderr="",
        )

    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/blastp" if name == "blastp" else None)
    monkeypatch.setattr(subprocess, "run", fake_run)

    search_rebase(
        query_faa=query,
        out_tsv=out_tsv,
        threads=1,
        mmseqs_db=tmp_path / "missing_mmseqs",
        blast_db=blast_prefix,
    )

    lines = out_tsv.read_text().splitlines()
    assert lines == [
        "query\ttarget\tevalue\tbits\tpident\tqcov\ttcov",
        "cand1\tM.EcoDcm\t0.0\t987\t100.000\t100\t100",
    ]


def _run_motif_calls(
    tmp_path: Path,
    *,
    search_rows: list[str],
    min_pident: float = 30.0,
    min_qcov: float = 50.0,
    min_tcov: float = 50.0,
    protein_map_rows: Optional[list[str]] = None,
    cluster_label_rows: Optional[list[str]] = None,
) -> tuple[Path, Path, Path]:
    candidates = tmp_path / "candidates.tsv"
    search = tmp_path / "search.tsv"
    enzymes = tmp_path / "enzymes.tsv"
    protein_map = tmp_path / "rebase_protein_map.tsv"
    cluster_labels = tmp_path / "rebase_cluster_labels.tsv"
    call_tsv = tmp_path / "cand1.tsv"
    consensus = tmp_path / "consensus.txt"
    pwm = tmp_path / "pwm.meme"

    _candidate_tsv(candidates)
    search.write_text("query\ttarget\tevalue\tbits\tpident\tqcov\ttcov\n" + "".join(search_rows))
    _enzymes_tsv(enzymes)
    if protein_map_rows is not None:
        protein_map.write_text(
            "protein_id\tenzyme_id\tmatched_name\tmotif_iupac\tmethylation\tmapping_method\n"
            + "".join(protein_map_rows)
        )
    if cluster_label_rows is not None:
        cluster_labels.write_text(
            "protein_id\tcluster_id\trepresentative_id\tcluster_size\tmotif_iupac\tmethylation\tconfidence\tlabel_method\tbest_subject_id\tbest_evalue\tbest_bits\tbest_pident\tbest_qcov\tbest_tcov\tsupport_count\n"
            + "".join(cluster_label_rows)
        )

    infer_motifs(
        candidates_tsv=candidates,
        search_tsv=search,
        enzymes_tsv=enzymes,
        out_calls=call_tsv,
        candidate_id="cand1",
        consensus_out=consensus,
        pwm_out=pwm,
        protein_map_tsv=protein_map if protein_map_rows is not None else None,
        cluster_labels_tsv=cluster_labels if cluster_label_rows is not None else None,
        max_evalue=1e-5,
        min_pident=min_pident,
        min_qcov=min_qcov,
        min_tcov=min_tcov,
        mmcif_cache_dir=tmp_path / "mmcif",
        mmcif_download=False,
        template_top_hits=3,
        template_min_alntmscore=0.5,
        panel_kmin=6,
        panel_kmax=8,
        panel_size=25,
    )
    return call_tsv, consensus, pwm


def _run_scan_and_qc(
    *,
    genome: Path,
    candidates: Path,
    call_tsv: Path,
    pwm: Path,
    fimo_tsv: Path,
    qc_json: Path,
    summary_tsv: Path,
) -> None:
    scan_and_summarize(
        genome_fa=genome,
        candidates_tsv=candidates,
        candidate_id="cand1",
        call_tsv=call_tsv,
        pwm_path=pwm,
        out_fimo=fimo_tsv,
        out_qc=qc_json,
        out_summary=summary_tsv,
        fimo_thresh=1e-4,
        bg_order=0,
    )


def test_motif_calls_candidate_mode_uses_rebase_cluster_transfer(tmp_path: Path) -> None:
    call_tsv, consensus, _pwm = _run_motif_calls(
        tmp_path,
        search_rows=["cand1\tM.OrfCluster1P\t0\t600\t100\t100\t100\n"],
        cluster_label_rows=[
            "M.OrfCluster1P\tcluster1\tM.OrfCluster1P\t18\tATGCAT\tm6A\thigh\tneighbor_duplicate_consensus\tM.AvaIII\t7.9e-108\t311\t55.435\t93\t93.5593\t1\n"
        ],
    )

    lines = call_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["hit_id"] == "M.AvaIII"
    assert row["motif_iupac"] == "ATGCAT"
    assert row["methylation"] == "m6A"
    assert row["mod_position"] == ""
    assert row["motif_class"] == "palindromic"
    assert row["motif_canonical_iupac"] == "ATGCAT"
    assert row["motif_reverse_complement_iupac"] == "ATGCAT"
    assert row["assignment_state"] == "linked"
    assert row["confidence"] == "high"
    assert row["method"] == "rebase_cluster_transfer"
    assert row["transfer_source_id"] == "M.OrfCluster1P"
    assert consensus.read_text() == "ATGCAT\n"


def test_motif_calls_candidate_mode_reports_no_motif(tmp_path: Path) -> None:
    call_tsv, consensus, pwm = _run_motif_calls(tmp_path, search_rows=[])

    lines = call_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["candidate_id"] == "cand1"
    assert row["motif_iupac"] == ""
    assert row["assignment_state"] == "unresolved"
    assert row["confidence"] == ""
    assert row["method"] == "no_motif"
    assert consensus.read_text() == "\n"
    assert "# No motif inferred for cand1" in pwm.read_text()


def test_motif_calls_candidate_mode_respects_foldseek_label_threshold(tmp_path: Path) -> None:
    candidates = tmp_path / "candidates.tsv"
    search = tmp_path / "search.tsv"
    enzymes = tmp_path / "enzymes.tsv"
    foldseek_hits = tmp_path / "foldseek_hits.tsv"
    foldseek_labels = tmp_path / "foldseek_labels.tsv"
    call_tsv = tmp_path / "cand1.tsv"
    consensus = tmp_path / "consensus.txt"
    pwm = tmp_path / "pwm.meme"

    _candidate_tsv(candidates)
    search.write_text("query\ttarget\tevalue\tbits\tpident\tqcov\ttcov\n")
    enzymes.write_text("enzyme_id\tnames\ttype\trecognition_seq_iupac\tmethylation\tevidence\n")
    foldseek_hits.write_text("query\ttarget\tevalue\tbits\talntmscore\ncand1\t1abc_A\t1e-10\t100\t0.6\n")
    foldseek_labels.write_text("target\tmotif_iupac\tmethylation\n1abc_A\tGATC\tm6A\n")

    infer_motifs(
        candidates_tsv=candidates,
        search_tsv=search,
        enzymes_tsv=enzymes,
        out_calls=call_tsv,
        candidate_id="cand1",
        consensus_out=consensus,
        pwm_out=pwm,
        foldseek_hits_tsv=foldseek_hits,
        foldseek_labels_tsv=foldseek_labels,
        max_evalue=1e-5,
        min_pident=30.0,
        min_qcov=50.0,
        min_tcov=50.0,
        mmcif_cache_dir=tmp_path / "mmcif",
        mmcif_download=False,
        template_top_hits=3,
        template_min_alntmscore=0.9,
        panel_kmin=6,
        panel_kmax=8,
        panel_size=25,
    )

    lines = call_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["motif_iupac"] == ""
    assert row["confidence"] == ""
    assert row["method"] == "no_motif"
    assert row["hit_id"] == ""
    assert consensus.read_text() == "\n"
    assert "# No motif inferred for cand1" in pwm.read_text()


def test_motif_calls_candidate_mode_uses_template_panel_without_labels(tmp_path: Path) -> None:
    candidates = tmp_path / "candidates.tsv"
    search = tmp_path / "search.tsv"
    enzymes = tmp_path / "enzymes.tsv"
    foldseek_hits = tmp_path / "foldseek_hits.tsv"
    mmcif_dir = tmp_path / "mmcif"
    call_tsv = tmp_path / "cand1.tsv"
    consensus = tmp_path / "consensus.txt"
    pwm = tmp_path / "pwm.meme"

    _candidate_tsv(candidates)
    search.write_text("query\ttarget\tevalue\tbits\tpident\tqcov\ttcov\n")
    enzymes.write_text("enzyme_id\tnames\ttype\trecognition_seq_iupac\tmethylation\tevidence\n")
    foldseek_hits.write_text("query\ttarget\tevalue\tbits\talntmscore\ncand1\t1abc_A\t1e-10\t100\t0.7\n")
    mmcif_dir.mkdir()
    (mmcif_dir / "1abc.cif").write_text(MMCIF_WITH_DNA)

    infer_motifs(
        candidates_tsv=candidates,
        search_tsv=search,
        enzymes_tsv=enzymes,
        out_calls=call_tsv,
        candidate_id="cand1",
        consensus_out=consensus,
        pwm_out=pwm,
        foldseek_hits_tsv=foldseek_hits,
        max_evalue=1e-5,
        min_pident=30.0,
        min_qcov=50.0,
        min_tcov=50.0,
        mmcif_mirror_dir=mmcif_dir,
        mmcif_cache_dir=tmp_path / "mmcif_cache",
        mmcif_download=False,
        template_top_hits=3,
        template_min_alntmscore=0.5,
        panel_kmin=4,
        panel_kmax=4,
        panel_size=1,
    )

    lines = call_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["motif_iupac"] == "GATC"
    assert row["method"] == "foldseek_template_panel"
    assert consensus.read_text() == "GATC\n"
    assert "MOTIF cand1" in pwm.read_text()


def test_motif_calls_candidate_mode_scales_fractional_mmseqs_metrics(tmp_path: Path) -> None:
    call_tsv, consensus, _pwm = _run_motif_calls(
        tmp_path,
        search_rows=["cand1\tEcoRI\t1e-40\t200\t0.55\t0.80\t0.75\n"],
        min_pident=50.0,
        min_qcov=50.0,
        min_tcov=50.0,
    )

    lines = call_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["hit_pident"] == "55"
    assert row["hit_qcov"] == "80"
    assert row["hit_tcov"] == "75"
    assert row["motif_iupac"] == "GAATTC"
    assert row["confidence"] == "medium"
    assert row["method"] == "rebase_homology"
    assert consensus.read_text() == "GAATTC\n"


def test_motif_calls_candidate_mode_uses_named_rebase_fallback_for_dcm(tmp_path: Path) -> None:
    call_tsv, consensus, _pwm = _run_motif_calls(
        tmp_path,
        search_rows=["cand1\tM.EcoPE144Dcm3P\t0\t987\t100\t100\t100\n"],
    )

    lines = call_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["motif_iupac"] == "CCWGG"
    assert row["methylation"] == "m5C"
    assert row["confidence"] == "high"
    assert row["method"] == "rebase_name_fallback"
    assert consensus.read_text() == "CCWGG\n"


def test_motif_calls_candidate_mode_uses_named_rebase_fallback_for_dam(tmp_path: Path) -> None:
    call_tsv, consensus, _pwm = _run_motif_calls(
        tmp_path,
        search_rows=["cand1\tM.EcoY204DamP\t0\t579\t99.64\t100\t100\n"],
    )

    lines = call_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["motif_iupac"] == "GATC"
    assert row["methylation"] == "m6A"
    assert row["confidence"] == "high"
    assert row["method"] == "rebase_name_fallback"
    assert consensus.read_text() == "GATC\n"


def test_motif_calls_candidate_mode_uses_rebase_protein_map(tmp_path: Path) -> None:
    call_tsv, consensus, _pwm = _run_motif_calls(
        tmp_path,
        search_rows=["cand1\tM.AatII\t1e-50\t300\t90\t100\t100\n"],
        protein_map_rows=["M.AatII\tAatII\tAatII\tGACGTC\tm6A\talias_strip_prefix\n"],
    )

    lines = call_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["motif_iupac"] == "GACGTC"
    assert row["methylation"] == "m6A"
    assert row["confidence"] == "high"
    assert row["method"] == "rebase_protein_map"
    assert consensus.read_text() == "GACGTC\n"


def test_motif_calls_candidate_mode_uses_rebase_consensus_when_best_hit_is_unmapped(tmp_path: Path) -> None:
    call_tsv, consensus, _pwm = _run_motif_calls(
        tmp_path,
        search_rows=[
            "cand1\tUnknownTop\t1e-80\t500\t95\t100\t100\n",
            "cand1\tM.AatII\t1e-60\t300\t90\t100\t100\n",
            "cand1\tM.BspRI\t1e-55\t280\t90\t100\t100\n",
        ],
        protein_map_rows=[
            "M.AatII\tAatII\tAatII\tGACGTC\tm6A\talias_strip_prefix\n",
            "M.BspRI\tBspRI\tBspRI\tGACGTC\tm6A\talias_strip_prefix\n",
        ],
    )

    lines = call_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["hit_id"] == "M.AatII"
    assert row["motif_iupac"] == "GACGTC"
    assert row["methylation"] == "m6A"
    assert row["confidence"] == "high"
    assert row["method"] == "rebase_consensus"
    assert consensus.read_text() == "GACGTC\n"


def test_motif_calls_candidate_mode_uses_rebase_family_transfer_for_orf_hits(
    tmp_path: Path,
    monkeypatch,
) -> None:
    candidates = tmp_path / "candidates.tsv"
    search = tmp_path / "search.tsv"
    enzymes = tmp_path / "enzymes.tsv"
    protein_map = tmp_path / "rebase_protein_map.tsv"
    rebase_proteins = tmp_path / "rebase_proteins.faa"
    rebase_motif_proteins = tmp_path / "rebase_motif_proteins.faa"
    call_tsv = tmp_path / "cand1.tsv"
    consensus = tmp_path / "consensus.txt"
    pwm = tmp_path / "pwm.meme"

    _candidate_tsv(candidates)
    search.write_text(
        "query\ttarget\tevalue\tbits\tpident\tqcov\ttcov\n"
        "cand1\tM.Orf1P\t0\t600\t100\t100\t100\n"
        "cand1\tM.Orf2P\t0\t600\t100\t100\t100\n"
    )
    _enzymes_tsv(enzymes)
    protein_map.write_text(
        "protein_id\tenzyme_id\tmatched_name\tmotif_iupac\tmethylation\tmapping_method\n"
        "M.AvaIII\tAvaIII\tAvaIII\tATGCAT\t5(6)\talias_strip_prefix\n"
    )
    rebase_proteins.write_text(
        ">M.Orf1P\nMPEPTIDE\n"
        ">M.Orf2P\nMPEPTIDE\n"
    )
    rebase_motif_proteins.write_text(
        ">M.AvaIII\nMPEPTIDE\n"
    )

    call_count = {"n": 0}

    def fake_run_counting(argv, capture_output, text):
        if argv[0] == "/usr/bin/blastp" and "-subject" in argv:
            call_count["n"] += 1
            query_path = Path(argv[argv.index("-query") + 1])
            query_text = query_path.read_text()
            query_id = query_text.splitlines()[0][1:]
            return subprocess.CompletedProcess(
                argv,
                0,
                stdout=f"{query_id}\tM.AvaIII\t1e-80\t250\t55\t95\t295\t280\n",
                stderr="",
            )
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/blastp" if name == "blastp" else None)
    monkeypatch.setattr(subprocess, "run", fake_run_counting)

    infer_motifs(
        candidates_tsv=candidates,
        search_tsv=search,
        enzymes_tsv=enzymes,
        out_calls=call_tsv,
        candidate_id="cand1",
        consensus_out=consensus,
        pwm_out=pwm,
        protein_map_tsv=protein_map,
        rebase_proteins=rebase_proteins,
        rebase_motif_proteins=rebase_motif_proteins,
        max_evalue=1e-5,
        min_pident=30.0,
        min_qcov=50.0,
        min_tcov=50.0,
        mmcif_cache_dir=tmp_path / "mmcif",
        mmcif_download=False,
        template_top_hits=3,
        template_min_alntmscore=0.5,
        panel_kmin=6,
        panel_kmax=8,
        panel_size=25,
    )

    lines = call_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["hit_id"] == "M.AvaIII"
    assert row["motif_iupac"] == "ATGCAT"
    assert row["methylation"] == "m6A"
    assert row["confidence"] == "medium"
    assert row["method"] == "rebase_family_transfer"
    assert row["transfer_source_id"] == "M.Orf1P"
    assert row["transfer_source_pident"] == "100"
    assert consensus.read_text() == "ATGCAT\n"
    assert call_count["n"] == 2


def test_motif_calls_candidate_mode_promotes_related_candidate_motif_for_unresolved_hit(
    tmp_path: Path,
    monkeypatch,
) -> None:
    candidates = tmp_path / "candidates.tsv"
    search = tmp_path / "search.tsv"
    enzymes = tmp_path / "enzymes.tsv"
    proteins = tmp_path / "mtase_candidates.faa"
    call_tsv = tmp_path / "cand1.tsv"
    consensus = tmp_path / "consensus.txt"
    pwm = tmp_path / "pwm.meme"

    candidates.write_text(
        "candidate_id\tcontig\tstart\tend\tstrand\tdomains\ttype_hint\tmethylation_hint\tneighborhood_features\tconfidence\trationale\n"
        "cand1\tctg1\t1\t100\t+\tpfam:PF00145(DNA_methylase):1e-80\ttype_II_or_orphan_like\tm5C\twindow_bp=10000\t1.0\tpfam_hit\n"
        "cand2\tctg1\t200\t300\t+\tpfam:PF00145(DNA_methylase):1e-80\ttype_II_or_orphan_like\tm5C\twindow_bp=10000\t1.0\tpfam_hit\n"
    )
    search.write_text(
        "query\ttarget\tevalue\tbits\tpident\tqcov\ttcov\n"
        "cand1\tM.SenDORF2P\t0\t994\t100\t100\t84.2756\n"
        "cand2\tM.EcoPE144Dcm3P\t0\t987\t100\t100\t100\n"
    )
    _enzymes_tsv(enzymes)
    proteins.write_text(">cand1\nMPEPTIDE\n>cand2\nMPEPTIDE\n")

    def fake_run(argv, capture_output, text):
        assert argv[0] == "/usr/bin/blastp"
        return subprocess.CompletedProcess(
            argv,
            0,
            stdout="cand1\tcand2\t0.0\t640\t68.707\t93\t477\t441\n",
            stderr="",
        )

    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/blastp" if name == "blastp" else None)
    monkeypatch.setattr(subprocess, "run", fake_run)

    infer_motifs(
        candidates_tsv=candidates,
        search_tsv=search,
        enzymes_tsv=enzymes,
        out_calls=call_tsv,
        candidate_id="cand1",
        consensus_out=consensus,
        pwm_out=pwm,
        candidate_proteins=proteins,
        max_evalue=1e-5,
        min_pident=30.0,
        min_qcov=50.0,
        min_tcov=50.0,
        mmcif_cache_dir=tmp_path / "mmcif",
        mmcif_download=False,
        template_top_hits=3,
        template_min_alntmscore=0.5,
        panel_kmin=6,
        panel_kmax=8,
        panel_size=25,
    )

    lines = call_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["motif_iupac"] == "CCWGG"
    assert row["methylation"] == "m5C"
    assert row["confidence"] == "medium"
    assert row["method"] == "candidate_homology_transfer"
    assert row["related_candidate_id"] == "cand2"
    assert row["related_motif_iupac"] == "CCWGG"
    assert row["related_methylation"] == "m5C"
    assert row["related_confidence"] == "medium"
    assert row["related_method"] == "candidate_homology_family"
    assert consensus.read_text() == "CCWGG\n"


def test_motif_calls_candidate_mode_keeps_related_candidate_annotation_only_when_methylation_conflicts(
    tmp_path: Path,
    monkeypatch,
) -> None:
    candidates = tmp_path / "candidates.tsv"
    search = tmp_path / "search.tsv"
    enzymes = tmp_path / "enzymes.tsv"
    proteins = tmp_path / "mtase_candidates.faa"
    call_tsv = tmp_path / "cand1.tsv"
    consensus = tmp_path / "consensus.txt"
    pwm = tmp_path / "pwm.meme"

    candidates.write_text(
        "candidate_id\tcontig\tstart\tend\tstrand\tdomains\ttype_hint\tmethylation_hint\tneighborhood_features\tconfidence\trationale\n"
        "cand1\tctg1\t1\t100\t+\tpfam:PF01555(N6_N4_Mtase):1e-80\ttype_II_or_orphan_like\tm6A/m4C\twindow_bp=10000\t1.0\tpfam_hit\n"
        "cand2\tctg1\t200\t300\t+\tpfam:PF00145(DNA_methylase):1e-80\ttype_II_or_orphan_like\tm5C\twindow_bp=10000\t1.0\tpfam_hit\n"
    )
    search.write_text(
        "query\ttarget\tevalue\tbits\tpident\tqcov\ttcov\n"
        "cand1\tM.SenDORF2P\t0\t994\t100\t100\t84.2756\n"
        "cand2\tM.EcoPE144Dcm3P\t0\t987\t100\t100\t100\n"
    )
    _enzymes_tsv(enzymes)
    proteins.write_text(">cand1\nMPEPTIDE\n>cand2\nMPEPTIDE\n")

    def fake_run(argv, capture_output, text):
        assert argv[0] == "/usr/bin/blastp"
        return subprocess.CompletedProcess(
            argv,
            0,
            stdout="cand1\tcand2\t0.0\t640\t68.707\t93\t477\t441\n",
            stderr="",
        )

    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/blastp" if name == "blastp" else None)
    monkeypatch.setattr(subprocess, "run", fake_run)

    infer_motifs(
        candidates_tsv=candidates,
        search_tsv=search,
        enzymes_tsv=enzymes,
        out_calls=call_tsv,
        candidate_id="cand1",
        consensus_out=consensus,
        pwm_out=pwm,
        candidate_proteins=proteins,
        max_evalue=1e-5,
        min_pident=30.0,
        min_qcov=50.0,
        min_tcov=50.0,
        mmcif_cache_dir=tmp_path / "mmcif",
        mmcif_download=False,
        template_top_hits=3,
        template_min_alntmscore=0.5,
        panel_kmin=6,
        panel_kmax=8,
        panel_size=25,
    )

    lines = call_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["motif_iupac"] == ""
    assert row["method"] == "rebase_hit_unmapped"
    assert row["related_candidate_id"] == "cand2"
    assert row["related_motif_iupac"] == "CCWGG"
    assert consensus.read_text() == "\n"


def test_motif_calls_candidate_mode_records_hint_for_weak_rebase_neighbor(
    tmp_path: Path,
    monkeypatch,
) -> None:
    candidates = tmp_path / "candidates.tsv"
    search = tmp_path / "search.tsv"
    enzymes = tmp_path / "enzymes.tsv"
    protein_map = tmp_path / "rebase_protein_map.tsv"
    rebase_proteins = tmp_path / "rebase_proteins.faa"
    rebase_motif_proteins = tmp_path / "rebase_motif_proteins.faa"
    call_tsv = tmp_path / "cand1.tsv"
    assignment_tsv = tmp_path / "assignments.tsv"
    consensus = tmp_path / "consensus.txt"
    pwm = tmp_path / "pwm.meme"

    _candidate_tsv(candidates)
    search.write_text(
        "query\ttarget\tevalue\tbits\tpident\tqcov\ttcov\n"
        "cand1\tM.SenWT5ORFDP\t3.64e-173\t478\t100\t100\t100\n"
    )
    _enzymes_tsv(enzymes)
    protein_map.write_text(
        "protein_id\tenzyme_id\tmatched_name\tmotif_iupac\tmethylation\tmapping_method\n"
        "M.NciI\tNciI\tNciI\tCCSGG\tm4C\talias_strip_prefix\n"
        "M.HinfI\tHinfI\tHinfI\tGANTC\tm6A\talias_strip_prefix\n"
    )
    rebase_proteins.write_text(">M.SenWT5ORFDP\nMPEPTIDE\n")
    rebase_motif_proteins.write_text(">M.NciI\nMPEPTIDE\n>M.HinfI\nMPEPTIDQ\n")

    def fake_run(argv, capture_output, text):
        if argv[0] == "/usr/bin/blastp" and "-subject" in argv:
            query_path = Path(argv[argv.index("-query") + 1])
            query_id = query_path.read_text().splitlines()[0][1:]
            return subprocess.CompletedProcess(
                argv,
                0,
                stdout=(
                    f"{query_id}\tM.NciI\t1.59e-18\t79.3\t25.738\t92\t462\t213\n"
                    f"{query_id}\tM.HinfI\t3.22e-17\t75.1\t24.231\t92\t424\t307\n"
                ),
                stderr="",
            )
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/blastp" if name == "blastp" else None)
    monkeypatch.setattr(subprocess, "run", fake_run)

    infer_motifs(
        candidates_tsv=candidates,
        search_tsv=search,
        enzymes_tsv=enzymes,
        out_calls=call_tsv,
        assignment_out=assignment_tsv,
        candidate_id="cand1",
        consensus_out=consensus,
        pwm_out=pwm,
        protein_map_tsv=protein_map,
        rebase_proteins=rebase_proteins,
        rebase_motif_proteins=rebase_motif_proteins,
        max_evalue=1e-5,
        min_pident=30.0,
        min_qcov=50.0,
        min_tcov=50.0,
        mmcif_cache_dir=tmp_path / "mmcif",
        mmcif_download=False,
        template_top_hits=3,
        template_min_alntmscore=0.5,
        panel_kmin=6,
        panel_kmax=8,
        panel_size=25,
    )

    lines = call_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["motif_iupac"] == ""
    assert row["assignment_state"] == "ambiguous"
    assert row["method"] == "rebase_hit_unmapped"
    assert row["hint_motif_iupac"] == "CCSGG"
    assert row["hint_methylation"] == "m4C"
    assert row["hint_confidence"] == "low"
    assert row["hint_method"] == "rebase_family_hint"
    assert row["hint_hit_id"] == "M.NciI"
    assert consensus.read_text() == "\n"

    assignment_lines = assignment_tsv.read_text().splitlines()
    assignment_header = assignment_lines[0].split("\t")
    assignment_row = dict(zip(assignment_header, assignment_lines[1].split("\t")))
    assert assignment_row["assignment_state"] == "ambiguous"
    assert assignment_row["assignment_role"] == "alternate"
    assert assignment_row["motif_iupac"] == "CCSGG"
    assert assignment_row["motif_canonical_iupac"] == "CCSGG"
    assert assignment_row["mod_position"] == ""
    assert assignment_row["method"] == "rebase_family_hint"


def test_scan_and_qc_candidate_mode_no_motif_writes_placeholder_outputs(
    tmp_path: Path,
    monkeypatch,
) -> None:
    genome = tmp_path / "genome.fna"
    candidates = tmp_path / "candidates.tsv"
    call_tsv = tmp_path / "cand1.tsv"
    pwm = tmp_path / "pwm.meme"
    fimo_tsv = tmp_path / "fimo" / "fimo.tsv"
    qc_json = tmp_path / "qc" / "qc.json"
    summary_tsv = tmp_path / "summary" / "cand1.tsv"

    genome.write_text(">ctg1\nACGTACGTACGTACGT\n")
    _candidate_tsv(candidates)
    call_tsv.write_text(
        "candidate_id\thit_id\thit_evalue\thit_bits\thit_pident\thit_qcov\thit_tcov\tmotif_iupac\tmethylation\tconfidence\tmethod\trelated_candidate_id\trelated_motif_iupac\trelated_methylation\trelated_confidence\trelated_method\trelated_pident\trelated_qcov\trelated_tcov\n"
        "cand1\t\t\t\t\t\t\t\t\t\tno_motif\tcand2\tCCWGG\tm5C\tmedium\tcandidate_homology_family\t68.707\t93\t92.4539\n"
    )
    pwm.write_text("MEME version 4\n")

    monkeypatch.setattr(
        shutil,
        "which",
        lambda name: "/usr/bin/fimo" if name == "fimo" else None,
    )

    _run_scan_and_qc(
        genome=genome,
        candidates=candidates,
        call_tsv=call_tsv,
        pwm=pwm,
        fimo_tsv=fimo_tsv,
        qc_json=qc_json,
        summary_tsv=summary_tsv,
    )

    assert fimo_tsv.read_text().startswith("motif_id\tmotif_alt_id\tsequence_name")
    qc = json.loads(qc_json.read_text())
    assert qc["warnings"] == ["no_motif"]
    lines = summary_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["candidate_methylation_hint"] == "m6A/m4C"
    assert row["candidate_confidence"] == "0.9"
    assert row["predicted_motif"] == ""
    assert row["predicted_motif_canonical"] == ""
    assert row["mod_position"] == ""
    assert row["motif_class"] == ""
    assert row["related_motif"] == "CCWGG"
    assert row["related_candidate_id"] == "cand2"
    assert row["related_method"] == "candidate_homology_family"
    assert row["assignment_state"] == "ambiguous"
    assert row["warnings"] == "no_motif"
    assert row["evidence_method"] == "no_motif"


def test_scan_and_qc_candidate_mode_uses_exact_iupac_fallback_when_fimo_returns_no_hits(
    tmp_path: Path,
    monkeypatch,
) -> None:
    genome = tmp_path / "genome.fna"
    candidates = tmp_path / "candidates.tsv"
    call_tsv = tmp_path / "cand1.tsv"
    pwm = tmp_path / "pwm.meme"
    fimo_tsv = tmp_path / "fimo" / "fimo.tsv"
    qc_json = tmp_path / "qc" / "qc.json"
    summary_tsv = tmp_path / "summary" / "cand1.tsv"

    genome.write_text(">ctg1\n" + ("A" * 2000) + "ACCAGG" + ("A" * 2000) + "CCTGG" + ("A" * 2000) + "\n")
    _candidate_tsv(candidates)
    call_tsv.write_text(
        "candidate_id\thit_id\thit_evalue\thit_bits\thit_pident\thit_qcov\thit_tcov\tmotif_iupac\tmethylation\tconfidence\tmethod\n"
        "cand1\tM.EcoDcm\t0\t987\t100\t100\t100\tCCWGG\tm5C\thigh\trebase_name_fallback\n"
    )
    pwm.write_text(
        "MEME version 4\n\n"
        "ALPHABET= ACGT\n\n"
        "strands: + -\n\n"
        "Background letter frequencies (from uniform background):\n"
        "A 0.25 C 0.25 G 0.25 T 0.25\n\n"
        "MOTIF cand1\n"
        "letter-probability matrix: alength= 4 w= 5 nsites= 20 E= 0\n"
        "0.000000 1.000000 0.000000 0.000000\n"
        "0.000000 1.000000 0.000000 0.000000\n"
        "0.500000 0.000000 0.000000 0.500000\n"
        "0.000000 0.000000 1.000000 0.000000\n"
        "0.000000 0.000000 1.000000 0.000000\n"
    )

    def fake_run(argv, capture_output, text):
        out_dir = Path(argv[argv.index("--oc") + 1])
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "fimo.tsv").write_text("# FIMO produced no rows at this threshold\n")
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(
        shutil,
        "which",
        lambda name: "/usr/bin/fimo" if name == "fimo" else None,
    )
    monkeypatch.setattr(subprocess, "run", fake_run)

    _run_scan_and_qc(
        genome=genome,
        candidates=candidates,
        call_tsv=call_tsv,
        pwm=pwm,
        fimo_tsv=fimo_tsv,
        qc_json=qc_json,
        summary_tsv=summary_tsv,
    )

    fimo_rows = [line for line in fimo_tsv.read_text().splitlines() if line and not line.startswith("#")]
    assert len(fimo_rows) == 5

    qc = json.loads(qc_json.read_text())
    assert qc["motif_iupac"] == "CCWGG"
    assert qc["num_sites"] == 2
    assert qc["warnings"] == []

    lines = summary_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["predicted_motif"] == "CCWGG"
    assert row["predicted_motif_canonical"] == "CCWGG"
    assert row["predicted_motif_reverse_complement"] == "CCWGG"
    assert row["mod_position"] == ""
    assert row["motif_class"] == "palindromic"
    assert row["assignment_state"] == "linked"
    assert row["num_sites"] == "2"
    assert row["warnings"] == ""


def test_scan_and_qc_candidate_mode_propagates_hint_fields_for_unresolved_candidate(
    tmp_path: Path,
    monkeypatch,
) -> None:
    genome = tmp_path / "genome.fna"
    candidates = tmp_path / "candidates.tsv"
    call_tsv = tmp_path / "cand1.tsv"
    pwm = tmp_path / "pwm.meme"
    fimo_tsv = tmp_path / "fimo" / "fimo.tsv"
    qc_json = tmp_path / "qc" / "qc.json"
    summary_tsv = tmp_path / "summary" / "cand1.tsv"

    genome.write_text(">ctg1\nACGTACGTACGTACGT\n")
    _candidate_tsv(candidates)
    call_tsv.write_text(
        "candidate_id\thit_id\thit_evalue\thit_bits\thit_pident\thit_qcov\thit_tcov\tmotif_iupac\tmethylation\tconfidence\tmethod\trelated_candidate_id\trelated_motif_iupac\trelated_methylation\trelated_confidence\trelated_method\trelated_pident\trelated_qcov\trelated_tcov\ttransfer_source_id\ttransfer_source_evalue\ttransfer_source_bits\ttransfer_source_pident\ttransfer_source_qcov\ttransfer_source_tcov\thint_motif_iupac\thint_methylation\thint_confidence\thint_method\thint_hit_id\thint_evalue\thint_bits\thint_pident\thint_qcov\thint_tcov\n"
        "cand1\tM.SenWT5ORFDP\t3.64e-173\t478\t100\t100\t100\t\t\t\trebase_hit_unmapped\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tCCSGG\tm4C\tlow\trebase_family_hint\tM.NciI\t1.59e-18\t79.3\t25.738\t92\t46.2\n"
    )
    pwm.write_text("MEME version 4\n")

    monkeypatch.setattr(
        shutil,
        "which",
        lambda name: "/usr/bin/fimo" if name == "fimo" else None,
    )

    _run_scan_and_qc(
        genome=genome,
        candidates=candidates,
        call_tsv=call_tsv,
        pwm=pwm,
        fimo_tsv=fimo_tsv,
        qc_json=qc_json,
        summary_tsv=summary_tsv,
    )

    lines = summary_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["predicted_motif"] == ""
    assert row["candidate_confidence"] == "0.9"
    assert row["hint_motif"] == "CCSGG"
    assert row["hint_methylation"] == "m4C"
    assert row["hint_hit_id"] == "M.NciI"
    assert row["hint_confidence"] == "low"
    assert row["hint_method"] == "rebase_family_hint"
    assert row["assignment_state"] == "ambiguous"
    assert row["warnings"] == "no_motif;hint_only"
    assert row["evidence_method"] == "rebase_hit_unmapped"


def test_scan_and_qc_candidate_mode_marks_low_support_unresolved_candidate(
    tmp_path: Path,
    monkeypatch,
) -> None:
    genome = tmp_path / "genome.fna"
    candidates = tmp_path / "candidates.tsv"
    call_tsv = tmp_path / "cand1.tsv"
    pwm = tmp_path / "pwm.meme"
    fimo_tsv = tmp_path / "fimo" / "fimo.tsv"
    qc_json = tmp_path / "qc" / "qc.json"
    summary_tsv = tmp_path / "summary" / "cand1.tsv"

    genome.write_text(">ctg1\nACGTACGTACGTACGT\n")
    candidates.write_text(
        "candidate_id\tcontig\tstart\tend\tstrand\tdomains\ttype_hint\tmethylation_hint\tneighborhood_features\tconfidence\trationale\n"
        "cand1\tctg1\t1\t100\t+\tpfam:PF01555(N6_N4_Mtase):4.1e-08\ttype_II_or_orphan_like\tm6A/m4C\twindow_bp=10000\t0.05\tpfam_hit\n"
    )
    call_tsv.write_text(
        "candidate_id\thit_id\thit_evalue\thit_bits\thit_pident\thit_qcov\thit_tcov\tmotif_iupac\tmethylation\tconfidence\tmethod\n"
        "cand1\t\t\t\t\t\t\t\t\t\tno_motif\n"
    )
    pwm.write_text("MEME version 4\n")

    monkeypatch.setattr(
        shutil,
        "which",
        lambda name: "/usr/bin/fimo" if name == "fimo" else None,
    )

    _run_scan_and_qc(
        genome=genome,
        candidates=candidates,
        call_tsv=call_tsv,
        pwm=pwm,
        fimo_tsv=fimo_tsv,
        qc_json=qc_json,
        summary_tsv=summary_tsv,
    )

    qc = json.loads(qc_json.read_text())
    assert qc["warnings"] == ["no_motif", "low_candidate_support"]

    lines = summary_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["candidate_confidence"] == "0.05"
    assert row["assignment_state"] == "unresolved"
    assert row["warnings"] == "no_motif;low_candidate_support"
    assert row["evidence_method"] == "no_motif"


def test_rebase_search_requires_blastp_for_blast_mode(tmp_path: Path, monkeypatch) -> None:
    query = tmp_path / "query.faa"
    out_tsv = tmp_path / "rebase.tsv"
    blast_prefix = tmp_path / "blast" / "rebase_proteins"

    query.write_text(">cand1\nMPEPTIDE\n")
    blast_prefix.parent.mkdir(parents=True)
    for ext in (".pin", ".psq", ".phr"):
        Path(str(blast_prefix) + ext).write_text("")

    monkeypatch.setattr(shutil, "which", lambda name: None if name == "blastp" else None)

    with pytest.raises(MtaseMotifError, match="blastp not found"):
        run_blast(
            query_faa=query,
            blast_db=blast_prefix,
            out_tsv=out_tsv,
            header="query\ttarget\tevalue\tbits\tpident\tqcov\ttcov\n",
            threads=1,
            max_target_seqs=10,
        )


def test_rebase_search_requires_mmseqs_for_mmseqs_mode(tmp_path: Path, monkeypatch) -> None:
    query = tmp_path / "query.faa"
    out_tsv = tmp_path / "rebase.tsv"

    query.write_text(">cand1\nMPEPTIDE\n")
    monkeypatch.setattr(shutil, "which", lambda name: None if name == "mmseqs" else None)

    with pytest.raises(MtaseMotifError, match="mmseqs not found"):
        run_mmseqs(
            query_faa=query,
            target_db=tmp_path / "mmseqs" / "rebase_proteins",
            out_tsv=out_tsv,
            header="query\ttarget\tevalue\tbits\tpident\tqcov\ttcov\n",
            threads=1,
        )


def test_scan_and_qc_rejects_malformed_candidates_tsv(
    tmp_path: Path,
    monkeypatch,
) -> None:
    genome = tmp_path / "genome.fna"
    candidates = tmp_path / "candidates.tsv"
    call_tsv = tmp_path / "cand1.tsv"
    pwm = tmp_path / "pwm.meme"
    fimo_tsv = tmp_path / "fimo" / "fimo.tsv"
    qc_json = tmp_path / "qc" / "qc.json"
    summary_tsv = tmp_path / "summary" / "cand1.tsv"

    genome.write_text(">ctg1\nACGTACGTACGTACGT\n")
    candidates.write_text("contig\tstart\tend\nctg1\t1\t100\n")
    call_tsv.write_text("candidate_id\tmotif_iupac\n")
    pwm.write_text("MEME version 4\n")

    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/fimo" if name == "fimo" else None)

    with pytest.raises(MtaseMotifError, match="Malformed candidates TSV"):
        scan_and_summarize(
            genome_fa=genome,
            candidates_tsv=candidates,
            candidate_id="cand1",
            call_tsv=call_tsv,
            pwm_path=pwm,
            out_fimo=fimo_tsv,
            out_qc=qc_json,
            out_summary=summary_tsv,
            fimo_thresh=1e-4,
            bg_order=0,
        )


def test_scan_and_qc_requires_candidate_mode_artifacts(
    tmp_path: Path,
    monkeypatch,
) -> None:
    genome = tmp_path / "genome.fna"
    candidates = tmp_path / "candidates.tsv"
    summary_tsv = tmp_path / "summary" / "cand1.tsv"

    genome.write_text(">ctg1\nACGTACGTACGTACGT\n")
    _candidate_tsv(candidates)
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/fimo" if name == "fimo" else None)

    with pytest.raises(MtaseMotifError, match="candidate_id mode requires"):
        scan_and_summarize(
            genome_fa=genome,
            candidates_tsv=candidates,
            candidate_id="cand1",
            call_tsv=None,
            pwm_path=None,
            out_fimo=None,
            out_qc=None,
            out_summary=summary_tsv,
            fimo_thresh=1e-4,
            bg_order=0,
        )


def test_run_foldseek_search_requires_foldseek_tool(
    tmp_path: Path,
    monkeypatch,
) -> None:
    structure_map = tmp_path / "structure_map.tsv"
    out_tsv = tmp_path / "foldseek_hits.tsv"
    structure_map.write_text("candidate_id\tstructure_path\ncand1\t/path/to/cand1.pdb\n")

    monkeypatch.setattr(shutil, "which", lambda name: None if name == "foldseek" else None)

    with pytest.raises(MtaseMotifError, match="foldseek not found"):
        run_foldseek_search(
            structure_map=structure_map,
            out_tsv=out_tsv,
            foldseek_db=tmp_path / "foldseek_db",
            threads=1,
        )


def test_find_structure_file_does_not_match_prefix_collisions(tmp_path: Path) -> None:
    structures_dir = tmp_path / "structures"
    structures_dir.mkdir()
    (structures_dir / "gene10_model.pdb").write_text("MODEL\n")
    assert find_structure_file(structures_dir, "gene1") is None

    (structures_dir / "gene1_model.pdb").write_text("MODEL\n")
    assert find_structure_file(structures_dir, "gene1") == structures_dir / "gene1_model.pdb"


def test_download_url_raises_on_empty_payload(tmp_path: Path, monkeypatch) -> None:
    dest = tmp_path / "downloaded.cif.gz"

    class EmptyResponse:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def read(self, size: int = -1) -> bytes:
            return b""

    monkeypatch.setattr("mtase_motif.structure_motif.template_panel.urllib.request.urlopen", lambda *args, **kwargs: EmptyResponse())

    with pytest.raises(MtaseMotifError, match="Downloaded empty file"):
        _download_url("https://example.org/template.cif.gz", dest)
