from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

import mtase_motif.cli as cli_module
from mtase_motif.cli import (
    FetchTarget,
    db_fetch,
    validate_pressed_hmm_db,
    validate_runtime_dependencies,
)
from mtase_motif.db.manifest import DbManifest, write_manifest
from mtase_motif.config import RunConfig
from mtase_motif.runner import build_run_plan, resolve_run_context
from mtase_motif.util import MtaseMotifError


def _seed_rebase_db(db_dir: Path) -> None:
    version_dir = db_dir / "rebase" / "20260317"
    version_dir.mkdir(parents=True)
    (version_dir / "rebase_enzymes.tsv").write_text(
        "enzyme_id\tnames\ttype\trecognition_seq_iupac\tmethylation\tevidence\n"
    )
    (version_dir / "rebase_proteins.faa").write_text(">M.Test\nMPEPTIDE\n")


def test_resolve_run_context_uses_default_structure_resources(tmp_path: Path) -> None:
    db_dir = tmp_path / "db"
    structures_dir = tmp_path / "structures"
    structures_dir.mkdir()
    _seed_rebase_db(db_dir)
    (db_dir / "hmms" / "pfam").mkdir(parents=True)
    (db_dir / "hmms" / "pfam" / "subset.hmm").write_text("HMMER3/f\n")
    (db_dir / "structures" / "pdb" / "mmcif").mkdir(parents=True)
    (db_dir / "structures" / "pdb" / "foldseek_db").mkdir(parents=True)
    (db_dir / "structures" / "pdb" / "labels.tsv").write_text("target\tmotif_iupac\n")

    cfg = RunConfig(
        genome_fasta=tmp_path / "genome.fa",
        db_dir=db_dir,
        out_dir=tmp_path / "out",
        structures_dir=structures_dir,
    )

    ctx = resolve_run_context(cfg)

    assert ctx.foldseek_db == db_dir / "structures" / "pdb" / "foldseek_db"
    assert ctx.foldseek_labels == db_dir / "structures" / "pdb" / "labels.tsv"
    assert ctx.mmcif_mirror_dir == db_dir / "structures" / "pdb" / "mmcif"
    assert ctx.foldseek_hits_tsv == tmp_path / "out" / "work" / "structures" / "foldseek_hits.tsv"


def test_build_run_plan_uses_native_steps(tmp_path: Path) -> None:
    db_dir = tmp_path / "db"
    _seed_rebase_db(db_dir)
    (db_dir / "hmms" / "pfam").mkdir(parents=True)
    (db_dir / "hmms" / "pfam" / "subset.hmm").write_text("HMMER3/f\n")

    ctx = resolve_run_context(
        RunConfig(
            genome_fasta=tmp_path / "genome.fa",
            db_dir=db_dir,
            out_dir=tmp_path / "out",
        )
    )

    plan = build_run_plan(ctx)

    assert plan[0] == "Call genes with Prodigal"
    assert "Search candidates against the REBASE protein index" in plan
    assert all("workflow engine" not in step.lower() for step in plan)


def test_validate_runtime_dependencies_skips_binary_checks_for_dry_run(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = RunConfig(
        genome_fasta=tmp_path / "genome.fa",
        db_dir=tmp_path / "db",
        out_dir=tmp_path / "out",
    )
    monkeypatch.setattr("mtase_motif.cli.which", lambda name: None)

    validate_runtime_dependencies(cfg, dry_run=True)


def test_validate_runtime_dependencies_requires_prodigal_for_real_run(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = RunConfig(
        genome_fasta=tmp_path / "genome.fa",
        db_dir=tmp_path / "db",
        out_dir=tmp_path / "out",
    )
    monkeypatch.setattr("mtase_motif.cli.which", lambda name: None)

    with pytest.raises(MtaseMotifError, match="prodigal not found"):
        validate_runtime_dependencies(cfg, dry_run=False)


def test_validate_runtime_dependencies_requires_foldseek_for_auto_discovered_db(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    db_dir = tmp_path / "db"
    structures_dir = tmp_path / "structures"
    structures_dir.mkdir()
    (db_dir / "structures" / "pdb" / "foldseek_db").mkdir(parents=True)

    cfg = RunConfig(
        genome_fasta=tmp_path / "genome.fa",
        db_dir=db_dir,
        out_dir=tmp_path / "out",
        structures_dir=structures_dir,
    )

    def fake_which(name: str):
        if name == "foldseek":
            return None
        return f"/usr/bin/{name}"

    monkeypatch.setattr("mtase_motif.cli.which", fake_which)

    with pytest.raises(MtaseMotifError, match="foldseek not found"):
        validate_runtime_dependencies(cfg, dry_run=False)


def test_validate_runtime_dependencies_requires_foldseek_db_for_labels_when_unavailable(
    tmp_path: Path,
) -> None:
    cfg = RunConfig(
        genome_fasta=tmp_path / "genome.fa",
        db_dir=tmp_path / "db",
        out_dir=tmp_path / "out",
        structures_dir=tmp_path / "structures",
        foldseek_labels=tmp_path / "labels.tsv",
    )

    with pytest.raises(MtaseMotifError, match="foldseek-labels requires a Foldseek DB"):
        validate_runtime_dependencies(cfg, dry_run=True)


def test_validate_runtime_dependencies_allows_labels_with_auto_discovered_foldseek_db(
    tmp_path: Path,
) -> None:
    db_dir = tmp_path / "db"
    structures_dir = tmp_path / "structures"
    structures_dir.mkdir()
    (db_dir / "structures" / "pdb" / "foldseek_db").mkdir(parents=True)

    cfg = RunConfig(
        genome_fasta=tmp_path / "genome.fa",
        db_dir=db_dir,
        out_dir=tmp_path / "out",
        structures_dir=structures_dir,
        foldseek_labels=tmp_path / "labels.tsv",
    )

    validate_runtime_dependencies(cfg, dry_run=True)


def test_validate_runtime_dependencies_requires_structures_dir_for_labels(
    tmp_path: Path,
) -> None:
    cfg = RunConfig(
        genome_fasta=tmp_path / "genome.fa",
        db_dir=tmp_path / "db",
        out_dir=tmp_path / "out",
        foldseek_labels=tmp_path / "labels.tsv",
    )

    with pytest.raises(MtaseMotifError, match="foldseek-labels requires --structures-dir"):
        validate_runtime_dependencies(cfg, dry_run=True)


def test_validate_pressed_hmm_db_requires_hmmpress_outputs(tmp_path: Path) -> None:
    hmm = tmp_path / "subset.hmm"
    hmm.write_text("HMMER3/f\n")

    with pytest.raises(MtaseMotifError, match="missing hmmpress index files"):
        validate_pressed_hmm_db(hmm, label="Pfam subset HMM")


def test_validate_pressed_hmm_db_accepts_pressed_hmm(tmp_path: Path) -> None:
    hmm = tmp_path / "subset.hmm"
    hmm.write_text("HMMER3/f\n")
    for ext in (".h3f", ".h3i", ".h3m", ".h3p"):
        (tmp_path / f"subset.hmm{ext}").write_text("")

    validate_pressed_hmm_db(hmm, label="Pfam subset HMM")


def test_run_config_normalized_rejects_invalid_panel_range(tmp_path: Path) -> None:
    cfg = RunConfig(
        genome_fasta=tmp_path / "genome.fa",
        db_dir=tmp_path / "db",
        out_dir=tmp_path / "out",
        panel_kmin=8,
        panel_kmax=6,
    )

    with pytest.raises(MtaseMotifError, match="panel-kmin"):
        cfg.normalized()


def test_main_wraps_unexpected_exceptions(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def boom() -> None:
        raise ValueError("boom")

    monkeypatch.setattr(cli_module, "app", boom)

    with pytest.raises(SystemExit) as excinfo:
        cli_module.main()

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert "Unexpected error" in captured.out
    assert "boom" in captured.out


def test_main_wraps_mtase_motif_errors(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def boom() -> None:
        raise MtaseMotifError("bad input")

    monkeypatch.setattr(cli_module, "app", boom)

    with pytest.raises(SystemExit) as excinfo:
        cli_module.main()

    assert excinfo.value.code == 2
    captured = capsys.readouterr()
    assert "Error:" in captured.out
    assert "bad input" in captured.out


def test_db_fetch_help_mentions_tigrfams_local_source_only() -> None:
    proc = subprocess.run(
        [sys.executable, "-m", "mtase_motif.cli", "db", "fetch", "--help"],
        capture_output=True,
        text=True,
    )

    assert proc.returncode == 0
    assert "local-source-only" in proc.stdout


def test_db_fetch_rejects_tigrfams_without_source(tmp_path: Path) -> None:
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    write_manifest(db_dir / "manifest.json", DbManifest())

    with pytest.raises(MtaseMotifError, match="requires --source"):
        db_fetch(FetchTarget.tigrfams, db_dir=db_dir, source=None)


def test_db_status_reports_provider_validation(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    write_manifest(db_dir / "manifest.json", DbManifest())

    cli_module.db_status(db_dir=db_dir)

    payload = json.loads(capsys.readouterr().out)
    assert payload["providers"]["pfam"]["valid"] is False
    assert payload["providers"]["tigrfams"]["valid"] is False
    assert payload["providers"]["rebase"]["valid"] is False
    assert "pfam not fetched" in payload["providers"]["pfam"]["validation_error"]
