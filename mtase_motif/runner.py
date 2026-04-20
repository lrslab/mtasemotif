from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from mtase_motif.candidate_discovery import build_candidates
from mtase_motif.config import RunConfig
from mtase_motif.motif_inference import infer_motifs
from mtase_motif.pipeline_support import (
    extract_candidate_proteins,
    run_foldseek_search,
    write_placeholder_gff,
    write_structure_map,
)
from mtase_motif.rebase_search import search_rebase
from mtase_motif.scanning import scan_and_summarize
from mtase_motif.util import MtaseMotifError, which


@dataclass(frozen=True)
class RunContext:
    cfg: RunConfig
    pfam_subset: Path
    tigr_subset: Optional[Path]
    rebase_version_dir: Path
    rebase_enzymes: Path
    rebase_proteins: Path
    rebase_protein_map: Optional[Path]
    rebase_cluster_labels: Optional[Path]
    rebase_motif_proteins: Optional[Path]
    mmseqs_db: Path
    blast_db: Path
    proteins_faa: Path
    gff_path: Path
    pfam_domtbl: Path
    pfam_log: Path
    tigr_domtbl: Path
    tigr_log: Path
    candidates_tsv: Path
    candidate_proteins_faa: Path
    rebase_search_tsv: Path
    motif_calls_tsv: Path
    motif_assignment_tsv: Path
    summary_tsv: Path
    structure_map_tsv: Optional[Path]
    foldseek_db: Optional[Path]
    foldseek_labels: Optional[Path]
    foldseek_hits_tsv: Optional[Path]
    mmcif_mirror_dir: Optional[Path]
    mmcif_cache_dir: Path


def resolve_run_context(cfg: RunConfig) -> RunContext:
    work_dir = cfg.out_dir / "work"
    pfam_subset = cfg.db_dir / "hmms" / "pfam" / "subset.hmm"
    tigr_path = cfg.db_dir / "hmms" / "tigrfams" / "subset.hmm"
    tigr_subset = tigr_path if tigr_path.exists() else None

    rebase_dir = cfg.db_dir / "rebase"
    rebase_versions = sorted([path for path in rebase_dir.iterdir() if path.is_dir()]) if rebase_dir.exists() else []
    if not rebase_versions:
        raise MtaseMotifError(
            f"No REBASE versions found under: {rebase_dir}\n"
            "Run: mtase-motif db fetch rebase"
        )
    rebase_version_dir = rebase_versions[-1]
    rebase_enzymes = rebase_version_dir / "rebase_enzymes.tsv"
    rebase_proteins = rebase_version_dir / "rebase_proteins.faa"
    rebase_protein_map = _optional_existing_path(rebase_version_dir / "rebase_protein_map.tsv")
    rebase_cluster_labels = _optional_existing_path(rebase_version_dir / "rebase_cluster_labels.tsv")
    rebase_motif_proteins = _optional_existing_path(rebase_version_dir / "rebase_motif_proteins.faa")
    mmseqs_db = rebase_version_dir / "mmseqs" / "rebase_proteins"
    blast_db = rebase_version_dir / "blast" / "rebase_proteins"

    if cfg.foldseek_db is not None:
        foldseek_db = cfg.foldseek_db
    else:
        default_foldseek_db = cfg.db_dir / "structures" / "pdb" / "foldseek_db"
        foldseek_db = default_foldseek_db if default_foldseek_db.exists() else None

    if cfg.foldseek_labels is not None:
        foldseek_labels = cfg.foldseek_labels
    elif foldseek_db is not None:
        default_foldseek_labels = foldseek_db.parent / "labels.tsv"
        foldseek_labels = default_foldseek_labels if default_foldseek_labels.exists() else None
    else:
        foldseek_labels = None

    if cfg.mmcif_mirror_dir is not None:
        mmcif_mirror_dir = cfg.mmcif_mirror_dir
    else:
        default_mmcif = cfg.db_dir / "structures" / "pdb" / "mmcif"
        mmcif_mirror_dir = default_mmcif if default_mmcif.exists() else None

    mmcif_cache_dir = cfg.mmcif_cache_dir or (work_dir / "structures" / "mmcif")
    proteins_faa = cfg.proteins_faa or (work_dir / "prodigal" / "proteins.faa")
    gff_path = work_dir / "prodigal" / "genes.gff"
    structure_map_tsv = (work_dir / "structures" / "structure_map.tsv") if cfg.structures_dir is not None else None
    foldseek_hits_tsv = (
        work_dir / "structures" / "foldseek_hits.tsv"
        if cfg.structures_dir is not None and foldseek_db is not None and foldseek_db.exists()
        else None
    )

    return RunContext(
        cfg=cfg,
        pfam_subset=pfam_subset,
        tigr_subset=tigr_subset,
        rebase_version_dir=rebase_version_dir,
        rebase_enzymes=rebase_enzymes,
        rebase_proteins=rebase_proteins,
        rebase_protein_map=rebase_protein_map,
        rebase_cluster_labels=rebase_cluster_labels,
        rebase_motif_proteins=rebase_motif_proteins,
        mmseqs_db=mmseqs_db,
        blast_db=blast_db,
        proteins_faa=proteins_faa,
        gff_path=gff_path,
        pfam_domtbl=work_dir / "hmmer" / "pfam.domtblout",
        pfam_log=work_dir / "hmmer" / "pfam.hmmscan.txt",
        tigr_domtbl=work_dir / "hmmer" / "tigrfams.domtblout",
        tigr_log=work_dir / "hmmer" / "tigrfams.hmmsearch.txt",
        candidates_tsv=cfg.out_dir / "mtase_candidates.tsv",
        candidate_proteins_faa=work_dir / "motifs" / "mtase_candidates.faa",
        rebase_search_tsv=work_dir / "motifs" / "rebase_search.tsv",
        motif_calls_tsv=cfg.out_dir / "motif_calls.tsv",
        motif_assignment_tsv=cfg.out_dir / "motif_assignment.tsv",
        summary_tsv=cfg.out_dir / "summary.tsv",
        structure_map_tsv=structure_map_tsv,
        foldseek_db=foldseek_db,
        foldseek_labels=foldseek_labels,
        foldseek_hits_tsv=foldseek_hits_tsv,
        mmcif_mirror_dir=mmcif_mirror_dir,
        mmcif_cache_dir=mmcif_cache_dir,
    )


def build_run_plan(ctx: RunContext) -> list[str]:
    steps = [
        "Use provided proteins" if ctx.cfg.proteins_faa is not None else "Call genes with Prodigal",
        "Scan proteins against the Pfam MTase subset",
        "Scan proteins against the TIGRFAMs subset" if ctx.tigr_subset is not None else "Skip TIGRFAMs scan",
        "Build and classify MTase candidates",
    ]
    if ctx.structure_map_tsv is not None:
        steps.append("Map candidate IDs to local structure files")
    if ctx.foldseek_hits_tsv is not None and ctx.foldseek_db is not None:
        steps.append(f"Search candidate structures against Foldseek DB {ctx.foldseek_db}")
    steps.extend(
        [
            "Extract candidate protein sequences",
            "Search candidates against the REBASE protein index",
            "Infer motifs from REBASE and optional structure evidence",
            "Scan the genome with FIMO and write QC + summary outputs",
        ]
    )
    return steps


def run_pipeline(
    ctx: RunContext,
    *,
    jobs: int = 1,
    dry_run: bool = False,
    logger: Optional[logging.Logger] = None,
) -> list[str]:
    plan = build_run_plan(ctx)
    if dry_run:
        return plan

    log = logger or logging.getLogger("mtase_motif.runner")
    _log_step(log, "gene_calling", ctx.cfg.proteins_faa or ctx.cfg.genome_fasta)
    if ctx.cfg.proteins_faa is None:
        _run_prodigal(genome=ctx.cfg.genome_fasta, proteins_out=ctx.proteins_faa, gff_out=ctx.gff_path)
    else:
        write_placeholder_gff(ctx.gff_path)

    _log_step(log, "hmmscan_pfam", ctx.pfam_subset)
    _run_hmmscan(
        hmm_path=ctx.pfam_subset,
        proteins_faa=ctx.proteins_faa,
        domtbl_out=ctx.pfam_domtbl,
        log_out=ctx.pfam_log,
        threads=jobs,
    )

    if ctx.tigr_subset is not None:
        _log_step(log, "hmmsearch_tigrfams", ctx.tigr_subset)
        _run_hmmsearch(
            hmm_path=ctx.tigr_subset,
            proteins_faa=ctx.proteins_faa,
            domtbl_out=ctx.tigr_domtbl,
            log_out=ctx.tigr_log,
            threads=jobs,
        )
    else:
        _write_skipped_tigr_outputs(ctx.tigr_domtbl, ctx.tigr_log)

    _log_step(log, "build_candidates", ctx.candidates_tsv)
    build_candidates(
        pfam_domtbl_path=ctx.pfam_domtbl,
        tigr_domtbl_path=ctx.tigr_domtbl,
        gff_path=ctx.gff_path,
        out_path=ctx.candidates_tsv,
        proteins_faa=ctx.proteins_faa,
        rebase_proteins=ctx.rebase_proteins,
    )

    if ctx.structure_map_tsv is not None and ctx.cfg.structures_dir is not None:
        _log_step(log, "structure_map", ctx.structure_map_tsv)
        write_structure_map(
            candidates_tsv=ctx.candidates_tsv,
            structures_dir=ctx.cfg.structures_dir,
            out_tsv=ctx.structure_map_tsv,
        )

    if ctx.foldseek_hits_tsv is not None and ctx.structure_map_tsv is not None and ctx.foldseek_db is not None:
        _log_step(log, "foldseek_search", ctx.foldseek_hits_tsv)
        run_foldseek_search(
            structure_map=ctx.structure_map_tsv,
            out_tsv=ctx.foldseek_hits_tsv,
            foldseek_db=ctx.foldseek_db,
            threads=jobs,
        )

    _log_step(log, "extract_candidate_proteins", ctx.candidate_proteins_faa)
    extract_candidate_proteins(
        proteins_faa=ctx.proteins_faa,
        candidates_tsv=ctx.candidates_tsv,
        out_faa=ctx.candidate_proteins_faa,
    )

    _log_step(log, "rebase_search", ctx.rebase_search_tsv)
    search_rebase(
        query_faa=ctx.candidate_proteins_faa,
        out_tsv=ctx.rebase_search_tsv,
        mmseqs_db=ctx.mmseqs_db,
        blast_db=ctx.blast_db,
        threads=jobs,
    )

    _log_step(log, "motif_inference", ctx.motif_calls_tsv)
    infer_motifs(
        candidates_tsv=ctx.candidates_tsv,
        search_tsv=ctx.rebase_search_tsv,
        enzymes_tsv=ctx.rebase_enzymes,
        out_calls=ctx.motif_calls_tsv,
        assignment_out=ctx.motif_assignment_tsv,
        candidate_proteins=ctx.candidate_proteins_faa,
        protein_map_tsv=ctx.rebase_protein_map,
        cluster_labels_tsv=ctx.rebase_cluster_labels,
        rebase_proteins=ctx.rebase_proteins,
        rebase_motif_proteins=ctx.rebase_motif_proteins,
        foldseek_hits_tsv=ctx.foldseek_hits_tsv,
        foldseek_labels_tsv=ctx.foldseek_labels,
        max_evalue=1e-5,
        min_pident=30.0,
        min_qcov=50.0,
        min_tcov=50.0,
        mmcif_mirror_dir=ctx.mmcif_mirror_dir,
        mmcif_cache_dir=ctx.mmcif_cache_dir,
        mmcif_download=ctx.cfg.mmcif_download,
        template_top_hits=ctx.cfg.template_top_hits,
        template_min_alntmscore=ctx.cfg.template_min_alntmscore,
        panel_kmin=ctx.cfg.panel_kmin,
        panel_kmax=ctx.cfg.panel_kmax,
        panel_size=ctx.cfg.panel_size,
    )

    _log_step(log, "scan_and_qc", ctx.summary_tsv)
    scan_and_summarize(
        genome_fa=ctx.cfg.genome_fasta,
        candidates_tsv=ctx.candidates_tsv,
        motifs_tsv=ctx.motif_calls_tsv,
        out_summary=ctx.summary_tsv,
        fimo_thresh=1e-4,
        bg_order=0,
    )
    return plan


def _log_step(logger: logging.Logger, step: str, target: Path) -> None:
    logger.info("pipeline_step=%s target=%s", step, target)


def _run_prodigal(*, genome: Path, proteins_out: Path, gff_out: Path) -> None:
    prodigal = which("prodigal")
    if prodigal is None:
        raise MtaseMotifError("prodigal not found in PATH")
    proteins_out.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        [
            prodigal,
            "-i",
            str(genome),
            "-a",
            str(proteins_out),
            "-o",
            str(gff_out),
            "-f",
            "gff",
            "-p",
            "single",
            "-q",
        ],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise MtaseMotifError(proc.stderr.strip() or "prodigal failed")


def _run_hmmscan(
    *,
    hmm_path: Path,
    proteins_faa: Path,
    domtbl_out: Path,
    log_out: Path,
    threads: int,
) -> None:
    hmmscan = which("hmmscan")
    if hmmscan is None:
        raise MtaseMotifError("hmmscan not found in PATH")
    domtbl_out.parent.mkdir(parents=True, exist_ok=True)
    with log_out.open("w") as log_handle:
        proc = subprocess.run(
            [
                hmmscan,
                "--noali",
                "--domtblout",
                str(domtbl_out),
                "--cpu",
                str(threads),
                str(hmm_path),
                str(proteins_faa),
            ],
            stdout=log_handle,
            stderr=subprocess.PIPE,
            text=True,
        )
    if proc.returncode != 0:
        raise MtaseMotifError(proc.stderr.strip() or "hmmscan failed")


def _run_hmmsearch(
    *,
    hmm_path: Path,
    proteins_faa: Path,
    domtbl_out: Path,
    log_out: Path,
    threads: int,
) -> None:
    hmmsearch = which("hmmsearch")
    if hmmsearch is None:
        raise MtaseMotifError("hmmsearch not found in PATH")
    domtbl_out.parent.mkdir(parents=True, exist_ok=True)
    with log_out.open("w") as log_handle:
        proc = subprocess.run(
            [
                hmmsearch,
                "--noali",
                "--domtblout",
                str(domtbl_out),
                "--cpu",
                str(threads),
                str(hmm_path),
                str(proteins_faa),
            ],
            stdout=log_handle,
            stderr=subprocess.PIPE,
            text=True,
        )
    if proc.returncode != 0:
        raise MtaseMotifError(proc.stderr.strip() or "hmmsearch failed")


def _write_skipped_tigr_outputs(domtbl_out: Path, log_out: Path) -> None:
    domtbl_out.parent.mkdir(parents=True, exist_ok=True)
    domtbl_out.write_text("# TIGRFAMs subset not found; skipped TIGRFAMs hmmsearch.\n")
    log_out.write_text("TIGRFAMs subset not found; skipped TIGRFAMs hmmsearch.\n")


def _optional_existing_path(path: Path) -> Optional[Path]:
    return path if path.exists() else None
