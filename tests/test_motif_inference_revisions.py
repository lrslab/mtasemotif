from __future__ import annotations

import subprocess
from pathlib import Path

from mtase_motif import motif_inference


def test_candidate_mode_only_builds_requested_candidate_when_related_transfer_is_unavailable(
    tmp_path: Path,
    monkeypatch,
) -> None:
    candidates = tmp_path / "candidates.tsv"
    search = tmp_path / "search.tsv"
    enzymes = tmp_path / "enzymes.tsv"
    calls = tmp_path / "calls.tsv"
    consensus = tmp_path / "consensus.txt"
    pwm = tmp_path / "pwm.meme"

    candidates.write_text(
        "candidate_id\tcontig\tstart\tend\tstrand\tdomains\ttype_hint\tmethylation_hint\tneighborhood_features\tconfidence\trationale\n"
        "cand1\tctg1\t1\t100\t+\tpfam\torphan_like\tm6A\twindow_bp=10000\t1.0\ttest\n"
        "cand2\tctg1\t200\t300\t+\tpfam\torphan_like\tm6A\twindow_bp=10000\t1.0\ttest\n"
    )
    search.write_text("query\ttarget\tevalue\tbits\tpident\tqcov\ttcov\n")
    enzymes.write_text("enzyme_id\trecognition_seq_iupac\tmethylation\n")

    seen: list[str] = []

    def fake_build_candidate_call(candidate_id: str, **_kwargs):
        seen.append(candidate_id)
        return (
            {
                "candidate_id": candidate_id,
                "motif_iupac": "",
                "methylation": "",
                "confidence": "",
                "method": "no_motif",
                "degraded_features": "",
            },
            None,
        )

    monkeypatch.setattr(motif_inference, "build_candidate_call", fake_build_candidate_call)
    monkeypatch.setattr(motif_inference, "annotate_related_candidate_calls", lambda *args, **kwargs: None)
    monkeypatch.setattr(motif_inference, "promote_related_candidate_calls", lambda *args, **kwargs: None)

    motif_inference.infer_motifs(
        candidates_tsv=candidates,
        search_tsv=search,
        enzymes_tsv=enzymes,
        out_calls=calls,
        candidate_id="cand1",
        consensus_out=consensus,
        pwm_out=pwm,
        mmcif_cache_dir=tmp_path / "mmcif",
    )

    assert seen == ["cand1"]


def test_related_candidate_transfer_batches_blast_queries(tmp_path: Path, monkeypatch) -> None:
    results_by_id = {
        "cand1": (
            {"candidate_id": "cand1", "motif_iupac": "", "methylation": "", "degraded_features": ""},
            None,
        ),
        "cand2": (
            {"candidate_id": "cand2", "motif_iupac": "", "methylation": "", "degraded_features": ""},
            None,
        ),
        "cand3": (
            {
                "candidate_id": "cand3",
                "motif_iupac": "CCWGG",
                "methylation": "m5C",
                "confidence": "high",
                "method": "rebase_homology",
                "degraded_features": "",
            },
            None,
        ),
    }
    candidate_sequences = {
        "cand1": "MPEPTIDE",
        "cand2": "MPEPTIDQ",
        "cand3": "MPEPTIDE",
    }
    calls: list[list[str]] = []

    def fake_run(argv, capture_output, text):
        calls.append(argv)
        return subprocess.CompletedProcess(
            argv,
            0,
            stdout=(
                "cand1\tcand3\t0.0\t600\t70\t90\t100\t100\n"
                "cand2\tcand3\t0.0\t610\t72\t91\t100\t100\n"
            ),
            stderr="",
        )

    monkeypatch.setattr(subprocess, "run", fake_run)

    motif_inference.annotate_related_candidate_calls(
        results_by_id,
        candidate_sequences,
        blastp_exe="/usr/bin/blastp",
    )

    assert len(calls) == 1
    assert results_by_id["cand1"][0]["related_candidate_id"] == "cand3"
    assert results_by_id["cand2"][0]["related_candidate_id"] == "cand3"


def test_candidate_mode_records_degraded_feature_when_related_transfer_is_skipped(
    tmp_path: Path,
    monkeypatch,
) -> None:
    candidates = tmp_path / "candidates.tsv"
    search = tmp_path / "search.tsv"
    enzymes = tmp_path / "enzymes.tsv"
    proteins = tmp_path / "proteins.faa"
    calls = tmp_path / "calls.tsv"
    consensus = tmp_path / "consensus.txt"
    pwm = tmp_path / "pwm.meme"

    candidates.write_text(
        "candidate_id\tcontig\tstart\tend\tstrand\tdomains\ttype_hint\tmethylation_hint\tneighborhood_features\tconfidence\trationale\n"
        "cand1\tctg1\t1\t100\t+\tpfam\torphan_like\tm5C\twindow_bp=10000\t1.0\ttest\n"
        "cand2\tctg1\t200\t300\t+\tpfam\torphan_like\tm5C\twindow_bp=10000\t1.0\ttest\n"
    )
    search.write_text("query\ttarget\tevalue\tbits\tpident\tqcov\ttcov\n")
    enzymes.write_text("enzyme_id\trecognition_seq_iupac\tmethylation\n")
    proteins.write_text(">cand1\nMPEPTIDE\n>cand2\nMPEPTIDQ\n")

    monkeypatch.setattr(motif_inference.shutil, "which", lambda name: None)

    motif_inference.infer_motifs(
        candidates_tsv=candidates,
        search_tsv=search,
        enzymes_tsv=enzymes,
        out_calls=calls,
        candidate_id="cand1",
        consensus_out=consensus,
        pwm_out=pwm,
        candidate_proteins=proteins,
        mmcif_cache_dir=tmp_path / "mmcif",
    )

    lines = calls.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert "candidate_homology_transfer_skipped_blastp_missing" in row["degraded_features"]
