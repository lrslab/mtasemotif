from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from mtase_motif.methylation import methylation_labels, normalize_methylation
from mtase_motif.motifs import (
    canonicalize_iupac_orientation,
    classify_iupac_motif,
    infer_mod_position,
    normalize_iupac,
    reverse_complement_iupac,
    write_meme_pwm,
)
from mtase_motif.structure_motif import foldseek as structure_foldseek
from mtase_motif.structure_motif import template_panel as structure_template_panel
from mtase_motif.util import MtaseMotifError


CALL_COLS = [
    "candidate_id",
    "hit_id",
    "hit_evalue",
    "hit_bits",
    "hit_pident",
    "hit_qcov",
    "hit_tcov",
    "motif_iupac",
    "methylation",
    "mod_position",
    "motif_class",
    "motif_canonical_iupac",
    "motif_reverse_complement_iupac",
    "assignment_state",
    "confidence",
    "method",
    "related_candidate_id",
    "related_motif_iupac",
    "related_methylation",
    "related_confidence",
    "related_method",
    "related_pident",
    "related_qcov",
    "related_tcov",
    "transfer_source_id",
    "transfer_source_evalue",
    "transfer_source_bits",
    "transfer_source_pident",
    "transfer_source_qcov",
    "transfer_source_tcov",
    "hint_motif_iupac",
    "hint_methylation",
    "hint_confidence",
    "hint_method",
    "hint_hit_id",
    "hint_evalue",
    "hint_bits",
    "hint_pident",
    "hint_qcov",
    "hint_tcov",
    "degraded_features",
]

ASSIGNMENT_COLS = [
    "candidate_id",
    "assignment_state",
    "assignment_role",
    "assignment_rank",
    "motif_iupac",
    "methylation",
    "mod_position",
    "motif_class",
    "motif_canonical_iupac",
    "motif_reverse_complement_iupac",
    "confidence",
    "method",
    "hit_id",
    "related_candidate_id",
    "transfer_source_id",
]


log = logging.getLogger(__name__)
_WARNED_FEATURES: set[str] = set()


@dataclass(frozen=True)
class Hit:
    target: str
    evalue: float
    bits: float
    pident: float
    qcov: float
    tcov: float


@dataclass(frozen=True)
class ProteinMapEntry:
    protein_id: str
    enzyme_id: str
    matched_name: str
    motif_iupac: str
    methylation: str
    mapping_method: str


@dataclass(frozen=True)
class RebaseClusterEntry:
    protein_id: str
    cluster_id: str
    representative_id: str
    cluster_size: int
    motif_iupac: str
    methylation: str
    confidence: str
    label_method: str
    best_subject_id: str
    best_evalue: Optional[float]
    best_bits: Optional[float]
    best_pident: Optional[float]
    best_qcov: Optional[float]
    best_tcov: Optional[float]
    support_count: int


@dataclass(frozen=True)
class ResolvedHit:
    hit: Hit
    motif_iupac: str
    methylation: str
    method: str


@dataclass(frozen=True)
class SubjectHit:
    subject_id: str
    evalue: float
    bits: float
    pident: float
    qcov: float
    tcov: float


@dataclass(frozen=True)
class RebaseFamilyRescue:
    source_hit: Hit
    mapped_subject_id: str
    motif_iupac: str
    methylation: str
    evalue: float
    bits: float
    pident: float
    qcov: float
    tcov: float


@dataclass(frozen=True)
class RebaseClusterRescue:
    source_hit: Hit
    entry: RebaseClusterEntry


@dataclass(frozen=True)
class RebaseHint:
    source_hit: Hit
    mapped_subject_id: str
    motif_iupac: str
    methylation: str
    evalue: float
    bits: float
    pident: float
    qcov: float
    tcov: float


@dataclass(frozen=True)
class CandidateFamilyHit:
    candidate_id: str
    motif_iupac: str
    methylation: str
    confidence: str
    method: str
    pident: float
    qcov: float
    tcov: float
    evalue: float
    bits: float


@dataclass(frozen=True)
class PairwiseHit:
    query_id: str
    subject_id: str
    evalue: float
    bits: float
    pident: float
    qcov: float
    tcov: float


def infer_motifs(
    *,
    candidates_tsv: Path,
    search_tsv: Path,
    enzymes_tsv: Path,
    out_calls: Path,
    assignment_out: Optional[Path] = None,
    candidate_id: str = "",
    consensus_out: Optional[Path] = None,
    pwm_out: Optional[Path] = None,
    candidate_proteins: Optional[Path] = None,
    protein_map_tsv: Optional[Path] = None,
    cluster_labels_tsv: Optional[Path] = None,
    rebase_proteins: Optional[Path] = None,
    rebase_motif_proteins: Optional[Path] = None,
    foldseek_hits_tsv: Optional[Path] = None,
    foldseek_labels_tsv: Optional[Path] = None,
    max_evalue: float = 1e-5,
    min_pident: float = 30.0,
    min_qcov: float = 50.0,
    min_tcov: float = 50.0,
    mmcif_mirror_dir: Optional[Path] = None,
    mmcif_cache_dir: Optional[Path] = None,
    mmcif_download: bool = False,
    template_top_hits: int = 3,
    template_min_alntmscore: float = 0.5,
    panel_kmin: int = 6,
    panel_kmax: int = 8,
    panel_size: int = 25,
) -> list[tuple[dict[str, str], Optional[structure_template_panel.TemplatePanelMotif]]]:
    candidate_mode = bool(candidate_id)
    if candidate_mode and (consensus_out is None or pwm_out is None):
        raise MtaseMotifError("candidate mode requires both consensus_out and pwm_out")
    out_calls.parent.mkdir(parents=True, exist_ok=True)
    if mmcif_cache_dir is None:
        mmcif_cache_dir = out_calls.parent / "work" / "structures" / "mmcif"

    candidate_rows = load_candidate_rows(candidates_tsv)
    all_candidate_ids = list(candidate_rows)
    if candidate_mode and candidate_id not in candidate_rows:
        raise MtaseMotifError(
            f"candidate_id {candidate_id!r} was not found in {candidates_tsv}"
        )
    target_candidate_ids = [candidate_id] if candidate_mode else all_candidate_ids
    blastp = shutil.which("blastp")
    related_transfer_possible = (
        candidate_proteins is not None and blastp is not None and len(all_candidate_ids) > 1
    )
    analysis_candidate_ids = all_candidate_ids if (not candidate_mode or related_transfer_possible) else target_candidate_ids

    hits_by_query = load_hits(search_tsv, selected_queries=set(analysis_candidate_ids))
    enzyme_to_motif = load_enzyme_motifs(enzymes_tsv)
    protein_to_motif = load_rebase_protein_map(protein_map_tsv)
    rebase_cluster_labels = load_rebase_cluster_labels(cluster_labels_tsv)
    rebase_target_ids = {
        hit.target
        for current_candidate_id in analysis_candidate_ids
        for hit in hits_by_query.get(current_candidate_id, [])
    }
    rebase_target_sequences = load_selected_fasta_sequences(rebase_proteins, rebase_target_ids)
    candidate_sequences = load_fasta_sequences(
        candidate_proteins,
        selected_ids=None if related_transfer_possible or not candidate_mode else set(target_candidate_ids),
    )
    if foldseek_hits_tsv is not None:
        foldseek_hits_by_query = structure_foldseek.load_foldseek_hits(foldseek_hits_tsv)
    else:
        foldseek_hits_by_query = {}
    if foldseek_labels_tsv is not None:
        foldseek_labels = structure_foldseek.load_foldseek_labels(foldseek_labels_tsv)
    else:
        foldseek_labels = {}

    results_by_id: dict[str, tuple[dict[str, str], Optional[structure_template_panel.TemplatePanelMotif]]] = {}
    for current_candidate_id in analysis_candidate_ids:
        row, panel = build_candidate_call(
            current_candidate_id,
            candidate_row=candidate_rows.get(current_candidate_id, {}),
            hits_by_query=hits_by_query,
            enzyme_to_motif=enzyme_to_motif,
            protein_to_motif=protein_to_motif,
            rebase_cluster_labels=rebase_cluster_labels,
            rebase_target_sequences=rebase_target_sequences,
            rebase_motif_proteins=rebase_motif_proteins,
            blastp_exe=blastp,
            foldseek_hits_by_query=foldseek_hits_by_query,
            foldseek_labels=foldseek_labels,
            max_evalue=max_evalue,
            min_pident=min_pident,
            min_qcov=min_qcov,
            min_tcov=min_tcov,
            mmcif_mirror_dir=mmcif_mirror_dir,
            mmcif_cache_dir=mmcif_cache_dir,
            mmcif_download=mmcif_download,
            template_top_hits=template_top_hits,
            template_min_alntmscore=template_min_alntmscore,
            panel_kmin=panel_kmin,
            panel_kmax=panel_kmax,
            panel_size=panel_size,
        )
        results_by_id[current_candidate_id] = (row, panel)

    annotate_related_candidate_calls(
        results_by_id,
        candidate_sequences,
        blastp_exe=blastp,
        target_ids=set(target_candidate_ids),
    )
    promote_related_candidate_calls(
        results_by_id,
        candidate_rows=candidate_rows,
        target_ids=set(target_candidate_ids),
    )

    if candidate_mode and candidate_proteins is not None and blastp is None and len(all_candidate_ids) > 1:
        row_and_panel = results_by_id.get(candidate_id)
        if row_and_panel is not None and not normalize_iupac(row_and_panel[0].get("motif_iupac", "")):
            _append_degraded_feature(
                row_and_panel[0], "candidate_homology_transfer_skipped_blastp_missing"
            )
            _log_warning_once(
                "candidate_related_blastp_missing",
                "blastp not found in PATH; skipping related-candidate motif transfer. "
                "Install BLAST+ to enable candidate-to-candidate rescue.",
            )

    results = [results_by_id[cid] for cid in target_candidate_ids if cid in results_by_id]
    for row, _panel in results:
        finalize_call_row(row)
    write_call_rows(out_calls, [row for row, _panel in results])
    if assignment_out is not None:
        write_assignment_rows(assignment_out, [row for row, _panel in results])

    if candidate_mode:
        if results:
            write_candidate_outputs(
                candidate_id,
                results[0][0],
                results[0][1],
                consensus_out,
                pwm_out,
            )
        else:
            write_empty_outputs(candidate_id, consensus_out, pwm_out)
    else:
        out_dir = out_calls.parent
        for row, panel in results:
            current_candidate_id = row["candidate_id"]
            candidate_dir = out_dir / current_candidate_id / "motif"
            write_candidate_outputs(
                current_candidate_id,
                row,
                panel,
                candidate_dir / "consensus.txt",
                candidate_dir / "pwm.meme",
            )
    return results


def build_candidate_call(
    candidate_id: str,
    *,
    candidate_row: dict[str, str],
    hits_by_query: dict[str, list[Hit]],
    enzyme_to_motif: dict[str, tuple[str, str]],
    protein_to_motif: dict[str, ProteinMapEntry],
    rebase_cluster_labels: dict[str, RebaseClusterEntry],
    rebase_target_sequences: dict[str, str],
    rebase_motif_proteins: Optional[Path],
    blastp_exe: Optional[str],
    foldseek_hits_by_query: dict[str, list[structure_foldseek.FoldseekHit]],
    foldseek_labels: dict[str, tuple[str, str]],
    max_evalue: float,
    min_pident: float,
    min_qcov: float,
    min_tcov: float,
    mmcif_mirror_dir: Optional[Path],
    mmcif_cache_dir: Path,
    mmcif_download: bool,
    template_top_hits: int,
    template_min_alntmscore: float,
    panel_kmin: int,
    panel_kmax: int,
    panel_size: int,
) -> tuple[dict[str, str], Optional[structure_template_panel.TemplatePanelMotif]]:
    method = "no_motif"
    confidence = ""
    panel = None
    degraded_features: list[str] = []

    passing_hits = choose_passing_hits(
        hits_by_query.get(candidate_id, []),
        max_evalue=max_evalue,
        min_pident=min_pident,
        min_qcov=min_qcov,
        min_tcov=min_tcov,
    )
    can_attempt_rebase_family_search = _can_attempt_rebase_family_search(
        passing_hits,
        enzyme_to_motif=enzyme_to_motif,
        protein_to_motif=protein_to_motif,
        rebase_target_sequences=rebase_target_sequences,
        rebase_motif_proteins=rebase_motif_proteins,
    )
    hit = passing_hits[0] if passing_hits else None
    hit_id = hit.target if hit else ""
    hit_evalue = f"{hit.evalue:g}" if hit else ""
    hit_bits = f"{hit.bits:g}" if hit else ""
    hit_pident = f"{hit.pident:g}" if hit else ""
    hit_qcov = f"{hit.qcov:g}" if hit else ""
    hit_tcov = f"{hit.tcov:g}" if hit else ""
    transfer_source_id = ""
    transfer_source_evalue = ""
    transfer_source_bits = ""
    transfer_source_pident = ""
    transfer_source_qcov = ""
    transfer_source_tcov = ""
    hint_motif_iupac = ""
    hint_methylation = ""
    hint_confidence = ""
    hint_method = ""
    hint_hit_id = ""
    hint_evalue = ""
    hint_bits = ""
    hint_pident = ""
    hint_qcov = ""
    hint_tcov = ""
    motif_iupac = ""
    methylation = ""

    resolved_hits: list[ResolvedHit] = []
    for passing_hit in passing_hits:
        resolved = resolve_hit_motif(
            passing_hit,
            enzyme_to_motif=enzyme_to_motif,
            protein_to_motif=protein_to_motif,
        )
        if resolved is not None:
            resolved_hits.append(resolved)

    if hit is not None:
        resolved_best = resolve_hit_motif(
            hit,
            enzyme_to_motif=enzyme_to_motif,
            protein_to_motif=protein_to_motif,
        )
        if resolved_best is not None:
            motif_iupac = resolved_best.motif_iupac
            methylation = normalize_methylation(resolved_best.methylation)
            confidence = confidence_label(hit)
            method = resolved_best.method
        else:
            consensus = choose_consensus_hit(resolved_hits)
            if consensus is not None:
                hit = consensus.hit
                hit_id = hit.target
                hit_evalue = f"{hit.evalue:g}"
                hit_bits = f"{hit.bits:g}"
                hit_pident = f"{hit.pident:g}"
                hit_qcov = f"{hit.qcov:g}"
                hit_tcov = f"{hit.tcov:g}"
                motif_iupac = consensus.motif_iupac
                methylation = normalize_methylation(consensus.methylation)
                confidence = confidence_label(hit)
                method = "rebase_consensus"
            else:
                cluster_rescue = infer_rebase_cluster_transfer(
                    passing_hits,
                    enzyme_to_motif=enzyme_to_motif,
                    protein_to_motif=protein_to_motif,
                    cluster_labels=rebase_cluster_labels,
                )
                if cluster_rescue is not None:
                    hit_id = cluster_rescue.entry.best_subject_id or cluster_rescue.entry.representative_id
                    hit_evalue = _format_optional_float(cluster_rescue.entry.best_evalue)
                    hit_bits = _format_optional_float(cluster_rescue.entry.best_bits)
                    hit_pident = _format_optional_float(cluster_rescue.entry.best_pident)
                    hit_qcov = _format_optional_float(cluster_rescue.entry.best_qcov)
                    hit_tcov = _format_optional_float(cluster_rescue.entry.best_tcov)
                    motif_iupac = cluster_rescue.entry.motif_iupac
                    methylation = normalize_methylation(cluster_rescue.entry.methylation)
                    confidence = cluster_rescue.entry.confidence
                    method = "rebase_cluster_transfer"
                    transfer_source_id = cluster_rescue.source_hit.target
                    transfer_source_evalue = f"{cluster_rescue.source_hit.evalue:g}"
                    transfer_source_bits = f"{cluster_rescue.source_hit.bits:g}"
                    transfer_source_pident = f"{cluster_rescue.source_hit.pident:g}"
                    transfer_source_qcov = f"{cluster_rescue.source_hit.qcov:g}"
                    transfer_source_tcov = f"{cluster_rescue.source_hit.tcov:g}"
                else:
                    if blastp_exe is None and can_attempt_rebase_family_search:
                        degraded_features.append("rebase_family_transfer_skipped_blastp_missing")
                        _log_warning_once(
                            "rebase_family_transfer_blastp_missing",
                            "blastp not found in PATH; skipping REBASE family transfer rescue. "
                            "Install BLAST+ to enable motif recovery from unmapped REBASE families.",
                        )
                    else:
                        rescue = infer_rebase_family_transfer(
                            passing_hits,
                            enzyme_to_motif=enzyme_to_motif,
                            protein_to_motif=protein_to_motif,
                            rebase_target_sequences=rebase_target_sequences,
                            rebase_motif_proteins=rebase_motif_proteins,
                            blastp_exe=blastp_exe,
                        )
                        if rescue is not None:
                            hit_id = rescue.mapped_subject_id
                            hit_evalue = f"{rescue.evalue:g}"
                            hit_bits = f"{rescue.bits:g}"
                            hit_pident = f"{rescue.pident:g}"
                            hit_qcov = f"{rescue.qcov:g}"
                            hit_tcov = f"{rescue.tcov:g}"
                            motif_iupac = rescue.motif_iupac
                            methylation = normalize_methylation(rescue.methylation)
                            confidence = confidence_from_rebase_family_rescue(rescue)
                            method = "rebase_family_transfer"
                            transfer_source_id = rescue.source_hit.target
                            transfer_source_evalue = f"{rescue.source_hit.evalue:g}"
                            transfer_source_bits = f"{rescue.source_hit.bits:g}"
                            transfer_source_pident = f"{rescue.source_hit.pident:g}"
                            transfer_source_qcov = f"{rescue.source_hit.qcov:g}"
                            transfer_source_tcov = f"{rescue.source_hit.tcov:g}"
                        else:
                            method = "rebase_hit_unmapped"

    if not motif_iupac and foldseek_hits_by_query and foldseek_labels:
        foldseek_hit = structure_foldseek.choose_best_foldseek_hit(
            foldseek_hits_by_query.get(candidate_id, []),
            min_alntmscore=template_min_alntmscore,
        )
        if foldseek_hit is not None:
            label = foldseek_labels.get(foldseek_hit.target)
            if label is not None and candidate_accepts_transferred_motif(candidate_row, label[1]):
                motif_iupac, methylation = label
                motif_iupac = normalize_iupac(motif_iupac)
                methylation = normalize_methylation(methylation)
                hit_id = foldseek_hit.target
                hit_evalue = f"{foldseek_hit.evalue:g}"
                hit_bits = f"{foldseek_hit.bits:g}"
                hit_pident = ""
                hit_qcov = ""
                hit_tcov = ""
                confidence = structure_foldseek.confidence_from_tmscore(foldseek_hit.alntmscore)
                method = "foldseek_structure"

    if not motif_iupac and foldseek_hits_by_query:
        panel = structure_template_panel.infer_template_panel_motif(
            foldseek_hits_by_query.get(candidate_id, []),
            mmcif_mirror_dir=mmcif_mirror_dir,
            mmcif_cache_dir=mmcif_cache_dir,
            allow_download=mmcif_download,
            top_hits=template_top_hits,
            min_alntmscore=template_min_alntmscore,
            kmin=panel_kmin,
            kmax=panel_kmax,
            panel_size=panel_size,
        )
        if panel is not None:
            motif_iupac = panel.motif_iupac
            methylation = ""
            hit_id = panel.best_target
            hit_evalue = f"{panel.best_evalue:g}"
            hit_bits = f"{panel.best_bits:g}"
            hit_pident = ""
            hit_qcov = ""
            hit_tcov = ""
            confidence = panel.confidence
            method = "foldseek_template_panel"

    if not motif_iupac:
        if blastp_exe is None and can_attempt_rebase_family_search:
            degraded_features.append("rebase_family_hint_skipped_blastp_missing")
            _log_warning_once(
                "rebase_family_hint_blastp_missing",
                "blastp not found in PATH; skipping REBASE family hint recovery. "
                "Install BLAST+ to enable weak-family motif hints.",
            )
        else:
            hint = infer_rebase_family_hint(
                passing_hits,
                enzyme_to_motif=enzyme_to_motif,
                protein_to_motif=protein_to_motif,
                rebase_target_sequences=rebase_target_sequences,
                rebase_motif_proteins=rebase_motif_proteins,
                blastp_exe=blastp_exe,
            )
            if hint is not None:
                hint_motif_iupac = hint.motif_iupac
                hint_methylation = normalize_methylation(hint.methylation)
                hint_confidence = confidence_from_rebase_hint(hint)
                hint_method = "rebase_family_hint"
                hint_hit_id = hint.mapped_subject_id
                hint_evalue = f"{hint.evalue:g}"
                hint_bits = f"{hint.bits:g}"
                hint_pident = f"{hint.pident:g}"
                hint_qcov = f"{hint.qcov:g}"
                hint_tcov = f"{hint.tcov:g}"

    row = {
        "candidate_id": candidate_id,
        "hit_id": hit_id,
        "hit_evalue": hit_evalue,
        "hit_bits": hit_bits,
        "hit_pident": hit_pident,
        "hit_qcov": hit_qcov,
        "hit_tcov": hit_tcov,
        "motif_iupac": motif_iupac,
        "methylation": methylation,
        "mod_position": "",
        "motif_class": "",
        "motif_canonical_iupac": "",
        "motif_reverse_complement_iupac": "",
        "assignment_state": "",
        "confidence": confidence if motif_iupac else "",
        "method": method,
        "related_candidate_id": "",
        "related_motif_iupac": "",
        "related_methylation": "",
        "related_confidence": "",
        "related_method": "",
        "related_pident": "",
        "related_qcov": "",
        "related_tcov": "",
        "transfer_source_id": transfer_source_id,
        "transfer_source_evalue": transfer_source_evalue,
        "transfer_source_bits": transfer_source_bits,
        "transfer_source_pident": transfer_source_pident,
        "transfer_source_qcov": transfer_source_qcov,
        "transfer_source_tcov": transfer_source_tcov,
        "hint_motif_iupac": hint_motif_iupac,
        "hint_methylation": hint_methylation,
        "hint_confidence": hint_confidence,
        "hint_method": hint_method,
        "hint_hit_id": hint_hit_id,
        "hint_evalue": hint_evalue,
        "hint_bits": hint_bits,
        "hint_pident": hint_pident,
        "hint_qcov": hint_qcov,
        "hint_tcov": hint_tcov,
        "degraded_features": ";".join(dict.fromkeys(degraded_features)),
    }
    return row, panel


def annotate_related_candidate_calls(
    results_by_id: dict[str, tuple[dict[str, str], Optional[structure_template_panel.TemplatePanelMotif]]],
    candidate_sequences: dict[str, str],
    *,
    blastp_exe: Optional[str] = None,
    target_ids: Optional[set[str]] = None,
) -> None:
    if not results_by_id or not candidate_sequences:
        return

    resolved_rows = {
        candidate_id: row
        for candidate_id, (row, _panel) in results_by_id.items()
        if normalize_iupac(row.get("motif_iupac", ""))
    }
    if not resolved_rows:
        return

    selected_ids = set(results_by_id) if target_ids is None else {cid for cid in target_ids if cid in results_by_id}
    unresolved_queries: dict[str, str] = {}
    for candidate_id in selected_ids:
        row, _panel = results_by_id[candidate_id]
        if normalize_iupac(row.get("motif_iupac", "")):
            continue
        query_seq = candidate_sequences.get(candidate_id, "")
        if not query_seq:
            continue
        unresolved_queries[candidate_id] = query_seq

    if not unresolved_queries:
        return

    subject_sequences = {
        resolved_id: candidate_sequences.get(resolved_id, "")
        for resolved_id in resolved_rows
        if candidate_sequences.get(resolved_id, "")
    }
    if not subject_sequences:
        return

    if blastp_exe is None:
        _log_warning_once(
            "candidate_related_blastp_missing",
            "blastp not found in PATH; skipping related-candidate motif transfer for %d unresolved candidate(s). "
            "Install BLAST+ to enable candidate-to-candidate rescue.",
            len(unresolved_queries),
        )
        for candidate_id in unresolved_queries:
            row, _panel = results_by_id[candidate_id]
            _append_degraded_feature(row, "candidate_homology_transfer_skipped_blastp_missing")
        return

    hits_by_query = run_related_candidate_blast(
        blastp_exe=blastp_exe,
        query_sequences=unresolved_queries,
        subject_sequences=subject_sequences,
    )
    for candidate_id, query_hits in hits_by_query.items():
        best = choose_best_related_candidate(
            candidate_id,
            related_hits=query_hits,
            resolved_rows=resolved_rows,
        )
        if best is None:
            continue
        row, _panel = results_by_id[candidate_id]
        row["related_candidate_id"] = best.candidate_id
        row["related_motif_iupac"] = best.motif_iupac
        row["related_methylation"] = best.methylation
        row["related_confidence"] = best.confidence
        row["related_method"] = best.method
        row["related_pident"] = f"{best.pident:g}"
        row["related_qcov"] = f"{best.qcov:g}"
        row["related_tcov"] = f"{best.tcov:g}"


def promote_related_candidate_calls(
    results_by_id: dict[str, tuple[dict[str, str], Optional[structure_template_panel.TemplatePanelMotif]]],
    *,
    candidate_rows: dict[str, dict[str, str]],
    target_ids: Optional[set[str]] = None,
) -> None:
    if not results_by_id:
        return

    selected_ids = set(results_by_id) if target_ids is None else {cid for cid in target_ids if cid in results_by_id}
    for candidate_id in selected_ids:
        row, _panel = results_by_id[candidate_id]
        if normalize_iupac(row.get("motif_iupac", "")):
            continue

        related_motif = normalize_iupac(row.get("related_motif_iupac", ""))
        if not related_motif:
            continue

        candidate_row = candidate_rows.get(candidate_id, {})
        related_methylation = row.get("related_methylation", "")
        if not candidate_accepts_transferred_motif(candidate_row, related_methylation):
            continue

        row["motif_iupac"] = related_motif
        row["methylation"] = normalize_methylation(related_methylation)
        row["confidence"] = row.get("related_confidence", "") or "low"
        row["method"] = "candidate_homology_transfer"


def choose_best_related_candidate(
    candidate_id: str,
    *,
    related_hits: list[PairwiseHit],
    resolved_rows: dict[str, dict[str, str]],
) -> Optional[CandidateFamilyHit]:
    best: Optional[CandidateFamilyHit] = None
    for pairwise in related_hits:
        resolved_id = pairwise.subject_id
        if resolved_id == candidate_id:
            continue
        resolved_row = resolved_rows.get(resolved_id)
        if resolved_row is None:
            continue
        if pairwise.pident < 60.0 or pairwise.qcov < 80.0 or pairwise.tcov < 80.0:
            continue
        related = CandidateFamilyHit(
            candidate_id=resolved_id,
            motif_iupac=normalize_iupac(resolved_row.get("motif_iupac", "")),
            methylation=resolved_row.get("methylation", ""),
            confidence=confidence_from_related_candidate(pairwise.pident, pairwise.qcov, pairwise.tcov),
            method="candidate_homology_family",
            pident=pairwise.pident,
            qcov=pairwise.qcov,
            tcov=pairwise.tcov,
            evalue=pairwise.evalue,
            bits=pairwise.bits,
        )
        if not related.motif_iupac:
            continue
        if best is None or (related.evalue, -related.bits, -related.pident, -related.qcov, -related.tcov) < (
            best.evalue,
            -best.bits,
            -best.pident,
            -best.qcov,
            -best.tcov,
        ):
            best = related
    return best


def run_related_candidate_blast(
    *,
    blastp_exe: str,
    query_sequences: dict[str, str],
    subject_sequences: dict[str, str],
) -> dict[str, list[PairwiseHit]]:
    if not query_sequences or not subject_sequences:
        return {}

    with tempfile.TemporaryDirectory(prefix="motif_related_") as tmpdir:
        tmpdir_path = Path(tmpdir)
        query_faa = tmpdir_path / "queries.faa"
        subject_faa = tmpdir_path / "subjects.faa"
        _write_fasta_records(query_faa, query_sequences)
        _write_fasta_records(subject_faa, subject_sequences)
        proc = subprocess.run(
            [
                blastp_exe,
                "-query",
                str(query_faa),
                "-subject",
                str(subject_faa),
                "-max_target_seqs",
                str(max(10, len(subject_sequences))),
                "-max_hsps",
                "1",
                "-outfmt",
                "6 qseqid sseqid evalue bitscore pident qcovs slen length",
            ],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            raise MtaseMotifError(
                proc.stderr.strip()
                or (
                    "blastp failed during related-candidate motif transfer "
                    f"for {query_faa} against {subject_faa}"
                )
            )

    hits_by_query: dict[str, list[PairwiseHit]] = {}
    for raw in proc.stdout.splitlines():
        if not raw.strip():
            continue
        parts = raw.split("\t")
        if len(parts) != 8:
            continue
        try:
            evalue = float(parts[2])
            bits = float(parts[3])
            pident = float(parts[4])
            qcov = float(parts[5])
            slen = float(parts[6])
            alen = float(parts[7])
        except ValueError:
            continue
        tcov = 0.0 if slen <= 0 else (100.0 * alen / slen)
        hit = PairwiseHit(
            query_id=parts[0],
            subject_id=parts[1],
            evalue=evalue,
            bits=bits,
            pident=pident,
            qcov=qcov,
            tcov=tcov,
        )
        hits_by_query.setdefault(hit.query_id, []).append(hit)
    return hits_by_query


def confidence_from_related_candidate(pident: float, qcov: float, tcov: float) -> str:
    if pident >= 70.0 and qcov >= 90.0 and tcov >= 90.0:
        return "high"
    if pident >= 60.0 and qcov >= 80.0 and tcov >= 80.0:
        return "medium"
    return "low"


def candidate_accepts_transferred_motif(candidate_row: dict[str, str], related_methylation: str) -> bool:
    candidate_labels = methylation_labels(candidate_row.get("methylation_hint", ""))
    related_labels = methylation_labels(related_methylation)
    if not candidate_labels or not related_labels:
        return True
    return bool(candidate_labels & related_labels)


def _append_degraded_feature(row: dict[str, str], feature: str) -> None:
    current = [item for item in row.get("degraded_features", "").split(";") if item]
    if feature not in current:
        current.append(feature)
    row["degraded_features"] = ";".join(current)


def _log_warning_once(key: str, message: str, *args: object) -> None:
    if key in _WARNED_FEATURES:
        return
    _WARNED_FEATURES.add(key)
    log.warning(message, *args)


def _write_fasta_records(path: Path, records: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for seq_id, sequence in records.items():
            handle.write(f">{seq_id}\n{sequence}\n")


def _validated_tsv_header(path: Path, header_line: str, *, required_cols: list[str]) -> list[str]:
    header = header_line.split("\t")
    missing = [col for col in required_cols if col not in header]
    if missing:
        raise MtaseMotifError(
            f"Malformed TSV at {path}: missing required columns {', '.join(missing)}"
        )
    return header


def infer_assignment_state(row: dict[str, str]) -> str:
    if normalize_iupac(row.get("motif_iupac", "")):
        return "linked"
    if normalize_iupac(row.get("related_motif_iupac", "")) or normalize_iupac(row.get("hint_motif_iupac", "")):
        return "ambiguous"
    return "unresolved"


def finalize_call_row(row: dict[str, str]) -> None:
    row["motif_iupac"] = normalize_iupac(row.get("motif_iupac", ""))
    row["methylation"] = normalize_methylation(row.get("methylation", ""))
    row["related_motif_iupac"] = normalize_iupac(row.get("related_motif_iupac", ""))
    row["related_methylation"] = normalize_methylation(row.get("related_methylation", ""))
    row["hint_motif_iupac"] = normalize_iupac(row.get("hint_motif_iupac", ""))
    row["hint_methylation"] = normalize_methylation(row.get("hint_methylation", ""))

    motif_iupac = row["motif_iupac"]
    if motif_iupac:
        row["motif_reverse_complement_iupac"] = reverse_complement_iupac(motif_iupac)
        row["motif_canonical_iupac"] = canonicalize_iupac_orientation(motif_iupac)
        row["motif_class"] = classify_iupac_motif(motif_iupac)
        mod_position = infer_mod_position(motif_iupac, row.get("methylation", ""))
        row["mod_position"] = str(mod_position) if mod_position is not None else ""
    else:
        row["motif_reverse_complement_iupac"] = ""
        row["motif_canonical_iupac"] = ""
        row["motif_class"] = ""
        row["mod_position"] = ""

    row["assignment_state"] = infer_assignment_state(row)


def _build_assignment_row(
    *,
    candidate_id: str,
    assignment_state: str,
    assignment_role: str,
    assignment_rank: int,
    motif_iupac: str,
    methylation: str,
    confidence: str,
    method: str,
    hit_id: str = "",
    related_candidate_id: str = "",
    transfer_source_id: str = "",
) -> dict[str, str]:
    motif_norm = normalize_iupac(motif_iupac)
    methylation_norm = normalize_methylation(methylation)
    mod_position = infer_mod_position(motif_norm, methylation_norm) if motif_norm else None
    return {
        "candidate_id": candidate_id,
        "assignment_state": assignment_state,
        "assignment_role": assignment_role,
        "assignment_rank": str(assignment_rank),
        "motif_iupac": motif_norm,
        "methylation": methylation_norm,
        "mod_position": str(mod_position) if mod_position is not None else "",
        "motif_class": classify_iupac_motif(motif_norm),
        "motif_canonical_iupac": canonicalize_iupac_orientation(motif_norm),
        "motif_reverse_complement_iupac": reverse_complement_iupac(motif_norm) if motif_norm else "",
        "confidence": confidence,
        "method": method,
        "hit_id": hit_id,
        "related_candidate_id": related_candidate_id,
        "transfer_source_id": transfer_source_id,
    }


def build_assignment_rows(row: dict[str, str]) -> list[dict[str, str]]:
    candidate_id = row.get("candidate_id", "")
    assignment_state = row.get("assignment_state", "") or infer_assignment_state(row)
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str, str, str, str, str]] = set()

    def add_assignment(
        *,
        assignment_role: str,
        motif_iupac: str,
        methylation: str,
        confidence: str,
        method: str,
        hit_id: str = "",
        related_candidate_id: str = "",
        transfer_source_id: str = "",
    ) -> None:
        motif_norm = normalize_iupac(motif_iupac)
        if not motif_norm:
            return
        key = (
            assignment_role,
            motif_norm,
            normalize_methylation(methylation),
            method,
            hit_id,
            related_candidate_id,
        )
        if key in seen:
            return
        seen.add(key)
        rows.append(
            _build_assignment_row(
                candidate_id=candidate_id,
                assignment_state=assignment_state,
                assignment_role=assignment_role,
                assignment_rank=len(rows) + 1,
                motif_iupac=motif_norm,
                methylation=methylation,
                confidence=confidence,
                method=method,
                hit_id=hit_id,
                related_candidate_id=related_candidate_id,
                transfer_source_id=transfer_source_id,
            )
        )

    primary_motif = normalize_iupac(row.get("motif_iupac", ""))
    primary_methylation = normalize_methylation(row.get("methylation", ""))
    if primary_motif:
        add_assignment(
            assignment_role="primary",
            motif_iupac=primary_motif,
            methylation=primary_methylation,
            confidence=row.get("confidence", ""),
            method=row.get("method", ""),
            hit_id=row.get("hit_id", ""),
            related_candidate_id=row.get("related_candidate_id", ""),
            transfer_source_id=row.get("transfer_source_id", ""),
        )

    related_motif = normalize_iupac(row.get("related_motif_iupac", ""))
    related_methylation = normalize_methylation(row.get("related_methylation", ""))
    if related_motif and (related_motif, related_methylation) != (primary_motif, primary_methylation):
        add_assignment(
            assignment_role="alternate",
            motif_iupac=related_motif,
            methylation=related_methylation,
            confidence=row.get("related_confidence", ""),
            method=row.get("related_method", ""),
            related_candidate_id=row.get("related_candidate_id", ""),
        )

    hint_motif = normalize_iupac(row.get("hint_motif_iupac", ""))
    hint_methylation = normalize_methylation(row.get("hint_methylation", ""))
    if hint_motif and (hint_motif, hint_methylation) != (primary_motif, primary_methylation):
        add_assignment(
            assignment_role="alternate",
            motif_iupac=hint_motif,
            methylation=hint_methylation,
            confidence=row.get("hint_confidence", ""),
            method=row.get("hint_method", ""),
            hit_id=row.get("hint_hit_id", ""),
        )

    if rows:
        return rows

    return [
        _build_assignment_row(
            candidate_id=candidate_id,
            assignment_state=assignment_state,
            assignment_role="primary",
            assignment_rank=1,
            motif_iupac="",
            methylation="",
            confidence="",
            method=row.get("method", ""),
            hit_id=row.get("hit_id", ""),
            related_candidate_id=row.get("related_candidate_id", ""),
            transfer_source_id=row.get("transfer_source_id", ""),
        )
    ]


def write_call_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("\t".join(CALL_COLS) + "\n")
        for row in rows:
            handle.write("\t".join(row.get(col, "") for col in CALL_COLS) + "\n")


def write_assignment_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    assignment_rows: list[dict[str, str]] = []
    for row in rows:
        assignment_rows.extend(build_assignment_rows(row))
    with path.open("w") as handle:
        handle.write("\t".join(ASSIGNMENT_COLS) + "\n")
        for row in assignment_rows:
            handle.write("\t".join(row.get(col, "") for col in ASSIGNMENT_COLS) + "\n")


def write_candidate_outputs(
    candidate_id: str,
    row: dict[str, str],
    panel: Optional[structure_template_panel.TemplatePanelMotif],
    consensus_path: Path,
    pwm_path: Path,
) -> None:
    motif_iupac = normalize_iupac(row.get("motif_iupac", ""))
    if panel is not None and motif_iupac:
        structure_template_panel.write_panel_outputs(consensus_path, pwm_path, candidate_id, panel)
        return
    if motif_iupac:
        consensus_path.parent.mkdir(parents=True, exist_ok=True)
        consensus_path.write_text(motif_iupac + "\n")
        write_meme_pwm(pwm_path, motif_id=candidate_id, motif_iupac=motif_iupac)
        return
    write_empty_outputs(candidate_id, consensus_path, pwm_path)


def write_empty_outputs(candidate_id: str, consensus_path: Path, pwm_path: Path) -> None:
    consensus_path.parent.mkdir(parents=True, exist_ok=True)
    consensus_path.write_text("\n")
    write_empty_meme(pwm_path, motif_id=candidate_id)


def write_empty_meme(path: Path, *, motif_id: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "MEME version 4\n\n"
        "ALPHABET= ACGT\n\n"
        "strands: + -\n\n"
        "Background letter frequencies (from uniform background):\n"
        "A 0.25 C 0.25 G 0.25 T 0.25\n\n"
        f"# No motif inferred for {motif_id}\n"
    )


def load_candidate_rows(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    lines = path.read_text().splitlines()
    if not lines:
        return rows

    header = _validated_tsv_header(path, lines[0], required_cols=["candidate_id"])
    for raw in lines[1:]:
        if not raw:
            continue
        parts = raw.split("\t")
        row = {header[i]: parts[i] if i < len(parts) else "" for i in range(len(header))}
        candidate_id = row.get("candidate_id", "")
        if candidate_id:
            rows[candidate_id] = row
    return rows


def load_fasta_sequences(
    path: Optional[Path],
    *,
    selected_ids: Optional[set[str]] = None,
) -> dict[str, str]:
    if path is None or not path.exists():
        return {}

    seqs: dict[str, str] = {}
    current_id: Optional[str] = None
    chunks: list[str] = []
    keep = selected_ids is None

    def flush() -> None:
        nonlocal current_id, chunks, keep
        if current_id is None:
            return
        if keep:
            seqs[current_id] = "".join(chunks).rstrip("*")
        current_id = None
        chunks = []
        keep = selected_ids is None

    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                current_id = line[1:].split()[0]
                keep = selected_ids is None or current_id in selected_ids
                continue
            if current_id is not None and keep:
                chunks.append(line)
    flush()
    return seqs


def load_selected_fasta_sequences(path: Optional[Path], selected_ids: set[str]) -> dict[str, str]:
    if path is None or not path.exists() or not selected_ids:
        return {}

    seqs: dict[str, str] = {}
    current_id: Optional[str] = None
    keep = False
    chunks: list[str] = []

    def flush() -> None:
        nonlocal current_id, keep, chunks
        if current_id is None:
            return
        if keep:
            seqs[current_id] = "".join(chunks).rstrip("*")
        current_id = None
        keep = False
        chunks = []

    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                current_id = line[1:].split()[0]
                keep = current_id in selected_ids
                continue
            if current_id is not None and keep:
                chunks.append(line)
    flush()
    return seqs


def load_hits(path: Path, *, selected_queries: Optional[set[str]] = None) -> dict[str, list[Hit]]:
    hits: dict[str, list[Hit]] = {}
    lines = path.read_text().splitlines()
    if not lines:
        return hits
    header = _validated_tsv_header(
        path,
        lines[0],
        required_cols=["query", "target", "evalue", "bits", "pident", "qcov", "tcov"],
    )
    query_idx = header.index("query")
    target_idx = header.index("target")
    evalue_idx = header.index("evalue")
    bits_idx = header.index("bits")
    pident_idx = header.index("pident")
    qcov_idx = header.index("qcov")
    tcov_idx = header.index("tcov")
    for raw in lines[1:]:
        if not raw.strip():
            continue
        parts = raw.split("\t")
        if len(parts) <= max(query_idx, target_idx, evalue_idx, bits_idx, pident_idx, qcov_idx, tcov_idx):
            continue
        query, target = parts[query_idx], parts[target_idx]
        if selected_queries is not None and query not in selected_queries:
            continue
        try:
            evalue = float(parts[evalue_idx])
            bits = float(parts[bits_idx])
            pident = _percent_if_fraction(float(parts[pident_idx]))
            qcov = _percent_if_fraction(float(parts[qcov_idx]))
            tcov = _percent_if_fraction(float(parts[tcov_idx]))
        except ValueError:
            continue
        hits.setdefault(query, []).append(
            Hit(target=target, evalue=evalue, bits=bits, pident=pident, qcov=qcov, tcov=tcov)
        )
    return hits


def load_enzyme_motifs(path: Path) -> dict[str, tuple[str, str]]:
    motif: dict[str, tuple[str, str]] = {}
    lines = path.read_text().splitlines()
    if not lines:
        return motif

    header = _validated_tsv_header(
        path,
        lines[0],
        required_cols=["enzyme_id", "recognition_seq_iupac", "methylation"],
    )
    for raw in lines[1:]:
        if not raw:
            continue
        parts = raw.split("\t")
        row = {header[i]: parts[i] if i < len(parts) else "" for i in range(len(header))}
        enzyme_id = row.get("enzyme_id", "")
        if not enzyme_id:
            continue
        motif[enzyme_id] = (
            row.get("recognition_seq_iupac", ""),
            normalize_methylation(row.get("methylation", "")),
        )
    return motif


def load_rebase_protein_map(path: Optional[Path]) -> dict[str, ProteinMapEntry]:
    mapping: dict[str, ProteinMapEntry] = {}
    if path is None or not path.exists():
        return mapping

    lines = path.read_text().splitlines()
    if not lines:
        return mapping
    header = _validated_tsv_header(
        path,
        lines[0],
        required_cols=["protein_id", "enzyme_id", "matched_name", "motif_iupac", "methylation", "mapping_method"],
    )
    for raw in lines[1:]:
        if not raw:
            continue
        parts = raw.split("\t")
        row = {header[i]: parts[i] if i < len(parts) else "" for i in range(len(header))}
        entry = ProteinMapEntry(
            protein_id=row.get("protein_id", ""),
            enzyme_id=row.get("enzyme_id", ""),
            matched_name=row.get("matched_name", ""),
            motif_iupac=normalize_iupac(row.get("motif_iupac", "")),
            methylation=normalize_methylation(row.get("methylation", "")),
            mapping_method=row.get("mapping_method", ""),
        )
        if entry.protein_id:
            mapping[entry.protein_id] = entry
    return mapping


def load_rebase_cluster_labels(path: Optional[Path]) -> dict[str, RebaseClusterEntry]:
    mapping: dict[str, RebaseClusterEntry] = {}
    if path is None or not path.exists():
        return mapping

    lines = path.read_text().splitlines()
    if not lines:
        return mapping
    header = _validated_tsv_header(
        path,
        lines[0],
        required_cols=[
            "protein_id",
            "cluster_id",
            "representative_id",
            "cluster_size",
            "motif_iupac",
            "methylation",
            "confidence",
            "label_method",
            "best_subject_id",
            "best_evalue",
            "best_bits",
            "best_pident",
            "best_qcov",
            "best_tcov",
            "support_count",
        ],
    )
    for raw in lines[1:]:
        if not raw:
            continue
        parts = raw.split("\t")
        row = {header[i]: parts[i] if i < len(parts) else "" for i in range(len(header))}
        protein_id = row.get("protein_id", "")
        if not protein_id:
            continue
        mapping[protein_id] = RebaseClusterEntry(
            protein_id=protein_id,
            cluster_id=row.get("cluster_id", ""),
            representative_id=row.get("representative_id", ""),
            cluster_size=_safe_int(row.get("cluster_size", "")),
            motif_iupac=normalize_iupac(row.get("motif_iupac", "")),
            methylation=normalize_methylation(row.get("methylation", "")),
            confidence=row.get("confidence", ""),
            label_method=row.get("label_method", ""),
            best_subject_id=row.get("best_subject_id", ""),
            best_evalue=_optional_float(row.get("best_evalue", "")),
            best_bits=_optional_float(row.get("best_bits", "")),
            best_pident=_optional_float(row.get("best_pident", "")),
            best_qcov=_optional_float(row.get("best_qcov", "")),
            best_tcov=_optional_float(row.get("best_tcov", "")),
            support_count=_safe_int(row.get("support_count", "")),
        )
    return mapping


def choose_passing_hits(
    hits: list[Hit],
    *,
    max_evalue: float,
    min_pident: float,
    min_qcov: float,
    min_tcov: float,
) -> list[Hit]:
    passing = [
        hit
        for hit in hits
        if hit.evalue <= max_evalue
        and hit.pident >= min_pident
        and hit.qcov >= min_qcov
        and hit.tcov >= min_tcov
    ]
    passing.sort(key=lambda hit: (hit.evalue, -hit.bits))
    return passing


def resolve_hit_motif(
    hit: Hit,
    *,
    enzyme_to_motif: dict[str, tuple[str, str]],
    protein_to_motif: dict[str, ProteinMapEntry],
) -> Optional[ResolvedHit]:
    motif_iupac, methylation = enzyme_to_motif.get(hit.target, ("", ""))
    motif_iupac = normalize_iupac(motif_iupac)
    if motif_iupac:
        return ResolvedHit(
            hit=hit,
            motif_iupac=motif_iupac,
            methylation=methylation,
            method="rebase_homology",
        )

    mapped = protein_to_motif.get(hit.target)
    if mapped is not None and mapped.motif_iupac:
        return ResolvedHit(
            hit=hit,
            motif_iupac=normalize_iupac(mapped.motif_iupac),
            methylation=mapped.methylation,
            method="rebase_protein_map",
        )

    fallback = infer_named_rebase_motif(hit.target)
    if fallback is None:
        return None
    return ResolvedHit(
        hit=hit,
        motif_iupac=normalize_iupac(fallback[0]),
        methylation=fallback[1],
        method="rebase_name_fallback",
    )


def choose_consensus_hit(resolved_hits: list[ResolvedHit]) -> Optional[ResolvedHit]:
    if len(resolved_hits) < 2:
        return None

    groups: dict[tuple[str, str], list[ResolvedHit]] = {}
    for resolved in resolved_hits:
        key = (resolved.motif_iupac, resolved.methylation)
        groups.setdefault(key, []).append(resolved)

    ranked = sorted(
        groups.items(),
        key=lambda item: (
            -len(item[1]),
            -sum(res.hit.bits for res in item[1]),
            item[0][0],
            item[0][1],
        ),
    )
    best_key, best_hits = ranked[0]
    if len(best_hits) < 2:
        return None

    representative = sorted(
        best_hits,
        key=lambda resolved: (
            resolved.hit.evalue,
            -resolved.hit.bits,
            -resolved.hit.pident,
            -resolved.hit.qcov,
            -resolved.hit.tcov,
        ),
    )[0]
    return ResolvedHit(
        hit=representative.hit,
        motif_iupac=best_key[0],
        methylation=best_key[1],
        method="rebase_consensus",
    )


def infer_rebase_cluster_transfer(
    passing_hits: list[Hit],
    *,
    enzyme_to_motif: dict[str, tuple[str, str]],
    protein_to_motif: dict[str, ProteinMapEntry],
    cluster_labels: dict[str, RebaseClusterEntry],
) -> Optional[RebaseClusterRescue]:
    rescues: list[RebaseClusterRescue] = []
    for hit in passing_hits:
        if resolve_hit_motif(hit, enzyme_to_motif=enzyme_to_motif, protein_to_motif=protein_to_motif) is not None:
            continue
        entry = cluster_labels.get(hit.target)
        if entry is None or not entry.motif_iupac:
            continue
        rescues.append(RebaseClusterRescue(source_hit=hit, entry=entry))
        if len(rescues) >= 5:
            break
    return choose_consensus_rebase_cluster_transfer(rescues)


def choose_consensus_rebase_cluster_transfer(
    rescues: list[RebaseClusterRescue],
) -> Optional[RebaseClusterRescue]:
    if not rescues:
        return None

    groups: dict[tuple[str, str], list[RebaseClusterRescue]] = {}
    for rescue in rescues:
        key = (rescue.entry.motif_iupac, rescue.entry.methylation)
        groups.setdefault(key, []).append(rescue)

    ranked = sorted(
        groups.items(),
        key=lambda item: (
            -len(item[1]),
            -max(max_cluster_confidence_score(rescue.entry.confidence) for rescue in item[1]),
            -sum(rescue.source_hit.bits for rescue in item[1]),
            item[0][0],
            item[0][1],
        ),
    )
    best_key, best_group = ranked[0]
    representative = sorted(
        best_group,
        key=lambda rescue: (
            rescue.source_hit.evalue,
            -rescue.source_hit.bits,
            -rescue.source_hit.pident,
            -rescue.source_hit.qcov,
            -rescue.source_hit.tcov,
        ),
    )[0]
    best_conf = max_cluster_confidence_score(representative.entry.confidence)
    if len(best_group) < 2 and best_conf < 2:
        return None

    return RebaseClusterRescue(
        source_hit=representative.source_hit,
        entry=RebaseClusterEntry(
            protein_id=representative.entry.protein_id,
            cluster_id=representative.entry.cluster_id,
            representative_id=representative.entry.representative_id,
            cluster_size=representative.entry.cluster_size,
            motif_iupac=best_key[0],
            methylation=best_key[1],
            confidence=representative.entry.confidence,
            label_method=representative.entry.label_method,
            best_subject_id=representative.entry.best_subject_id,
            best_evalue=representative.entry.best_evalue,
            best_bits=representative.entry.best_bits,
            best_pident=representative.entry.best_pident,
            best_qcov=representative.entry.best_qcov,
            best_tcov=representative.entry.best_tcov,
            support_count=representative.entry.support_count,
        ),
    )


def infer_rebase_family_transfer(
    passing_hits: list[Hit],
    *,
    enzyme_to_motif: dict[str, tuple[str, str]],
    protein_to_motif: dict[str, ProteinMapEntry],
    rebase_target_sequences: dict[str, str],
    rebase_motif_proteins: Optional[Path],
    blastp_exe: Optional[str],
) -> Optional[RebaseFamilyRescue]:
    if not passing_hits or rebase_motif_proteins is None or not rebase_motif_proteins.exists():
        return None

    if blastp_exe is None:
        return None

    rescues: list[RebaseFamilyRescue] = []
    unresolved_seen = 0
    for hit in passing_hits:
        if resolve_hit_motif(hit, enzyme_to_motif=enzyme_to_motif, protein_to_motif=protein_to_motif) is not None:
            continue

        query_seq = rebase_target_sequences.get(hit.target, "")
        if not query_seq:
            continue

        rescue = choose_best_rebase_family_neighbor(
            blastp_exe=blastp_exe,
            source_hit=hit,
            query_seq=query_seq,
            subject_faa=rebase_motif_proteins,
            protein_to_motif=protein_to_motif,
        )
        if rescue is not None:
            rescues.append(rescue)

        unresolved_seen += 1
        if unresolved_seen >= 5:
            break

    return choose_consensus_rebase_family_rescue(rescues)


def infer_rebase_family_hint(
    passing_hits: list[Hit],
    *,
    enzyme_to_motif: dict[str, tuple[str, str]],
    protein_to_motif: dict[str, ProteinMapEntry],
    rebase_target_sequences: dict[str, str],
    rebase_motif_proteins: Optional[Path],
    blastp_exe: Optional[str],
) -> Optional[RebaseHint]:
    if not passing_hits or rebase_motif_proteins is None or not rebase_motif_proteins.exists():
        return None

    if blastp_exe is None:
        return None

    hints: list[RebaseHint] = []
    unresolved_seen = 0
    for hit in passing_hits:
        if resolve_hit_motif(hit, enzyme_to_motif=enzyme_to_motif, protein_to_motif=protein_to_motif) is not None:
            continue

        query_seq = rebase_target_sequences.get(hit.target, "")
        if not query_seq:
            continue

        hint = choose_best_rebase_hint_neighbor(
            blastp_exe=blastp_exe,
            source_hit=hit,
            query_seq=query_seq,
            subject_faa=rebase_motif_proteins,
            protein_to_motif=protein_to_motif,
        )
        if hint is not None:
            hints.append(hint)

        unresolved_seen += 1
        if unresolved_seen >= 5:
            break

    return choose_consensus_rebase_hint(hints)


def choose_best_rebase_family_neighbor(
    *,
    blastp_exe: str,
    source_hit: Hit,
    query_seq: str,
    subject_faa: Path,
    protein_to_motif: dict[str, ProteinMapEntry],
) -> Optional[RebaseFamilyRescue]:
    if not query_seq:
        return None

    subject_hits = run_blast_to_subject_fasta(
        blastp_exe=blastp_exe,
        query_id=source_hit.target,
        query_seq=query_seq,
        subject_faa=subject_faa,
        max_target_seqs=10,
    )
    best: Optional[RebaseFamilyRescue] = None
    for subject_hit in subject_hits:
        if (
            subject_hit.evalue > 1e-40
            or subject_hit.pident < 50.0
            or subject_hit.qcov < 80.0
            or subject_hit.tcov < 80.0
        ):
            continue

        mapped = protein_to_motif.get(subject_hit.subject_id)
        if mapped is None or not mapped.motif_iupac:
            continue

        rescue = RebaseFamilyRescue(
            source_hit=source_hit,
            mapped_subject_id=subject_hit.subject_id,
            motif_iupac=normalize_iupac(mapped.motif_iupac),
            methylation=mapped.methylation,
            evalue=subject_hit.evalue,
            bits=subject_hit.bits,
            pident=subject_hit.pident,
            qcov=subject_hit.qcov,
            tcov=subject_hit.tcov,
        )
        if best is None or (
            rescue.evalue,
            -rescue.bits,
            -rescue.pident,
            -rescue.qcov,
            -rescue.tcov,
        ) < (
            best.evalue,
            -best.bits,
            -best.pident,
            -best.qcov,
            -best.tcov,
        ):
            best = rescue
    return best


def choose_best_rebase_hint_neighbor(
    *,
    blastp_exe: str,
    source_hit: Hit,
    query_seq: str,
    subject_faa: Path,
    protein_to_motif: dict[str, ProteinMapEntry],
) -> Optional[RebaseHint]:
    if not query_seq:
        return None

    subject_hits = run_blast_to_subject_fasta(
        blastp_exe=blastp_exe,
        query_id=source_hit.target,
        query_seq=query_seq,
        subject_faa=subject_faa,
        max_target_seqs=10,
    )
    candidates: list[RebaseHint] = []
    for subject_hit in subject_hits:
        if (
            subject_hit.evalue > 1e-12
            or subject_hit.bits < 70.0
            or subject_hit.pident < 20.0
            or subject_hit.qcov < 85.0
            or subject_hit.tcov < 40.0
        ):
            continue

        mapped = protein_to_motif.get(subject_hit.subject_id)
        if mapped is None or not mapped.motif_iupac:
            continue

        candidates.append(
            RebaseHint(
                source_hit=source_hit,
                mapped_subject_id=subject_hit.subject_id,
                motif_iupac=normalize_iupac(mapped.motif_iupac),
                methylation=mapped.methylation,
                evalue=subject_hit.evalue,
                bits=subject_hit.bits,
                pident=subject_hit.pident,
                qcov=subject_hit.qcov,
                tcov=subject_hit.tcov,
            )
        )

    if not candidates:
        return None

    candidates.sort(
        key=lambda hint: (
            hint.evalue,
            -hint.bits,
            -hint.pident,
            -hint.qcov,
            -hint.tcov,
        )
    )
    best = candidates[0]
    second = candidates[1] if len(candidates) > 1 else None
    if second is not None:
        same_label = (
            best.motif_iupac == second.motif_iupac
            and normalize_methylation(best.methylation) == normalize_methylation(second.methylation)
        )
        if not same_label and best.bits < second.bits * 1.03:
            return None
    return best


def run_blast_to_subject_fasta(
    *,
    blastp_exe: str,
    query_id: str,
    query_seq: str,
    subject_faa: Path,
    max_target_seqs: int,
) -> list[SubjectHit]:
    with tempfile.TemporaryDirectory(prefix="motif_rebase_family_") as tmpdir:
        query_faa = Path(tmpdir) / "query.faa"
        query_faa.write_text(f">{query_id}\n{query_seq}\n")
        proc = subprocess.run(
            [
                blastp_exe,
                "-query",
                str(query_faa),
                "-subject",
                str(subject_faa),
                "-max_target_seqs",
                str(max_target_seqs),
                "-max_hsps",
                "1",
                "-outfmt",
                "6 qseqid sseqid evalue bitscore pident qcovs slen length",
            ],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            raise MtaseMotifError(
                proc.stderr.strip()
                or (
                    "blastp failed during REBASE family rescue "
                    f"for {query_faa} against {subject_faa}"
                )
            )

    hits: list[SubjectHit] = []
    for raw in proc.stdout.splitlines():
        if not raw.strip():
            continue
        parts = raw.split("\t")
        if len(parts) != 8:
            continue
        try:
            evalue = float(parts[2])
            bits = float(parts[3])
            pident = float(parts[4])
            qcov = float(parts[5])
            slen = float(parts[6])
            alen = float(parts[7])
        except ValueError:
            continue
        tcov = 0.0 if slen <= 0 else (100.0 * alen / slen)
        hits.append(
            SubjectHit(
                subject_id=parts[1],
                evalue=evalue,
                bits=bits,
                pident=pident,
                qcov=qcov,
                tcov=tcov,
            )
        )
    return hits


def _can_attempt_rebase_family_search(
    passing_hits: list[Hit],
    *,
    enzyme_to_motif: dict[str, tuple[str, str]],
    protein_to_motif: dict[str, ProteinMapEntry],
    rebase_target_sequences: dict[str, str],
    rebase_motif_proteins: Optional[Path],
) -> bool:
    if rebase_motif_proteins is None or not rebase_motif_proteins.exists():
        return False
    unresolved_seen = 0
    for hit in passing_hits:
        if resolve_hit_motif(hit, enzyme_to_motif=enzyme_to_motif, protein_to_motif=protein_to_motif) is not None:
            continue
        if rebase_target_sequences.get(hit.target, ""):
            return True
        unresolved_seen += 1
        if unresolved_seen >= 5:
            break
    return False


def choose_consensus_rebase_family_rescue(
    rescues: list[RebaseFamilyRescue],
) -> Optional[RebaseFamilyRescue]:
    if len(rescues) < 2:
        return None

    groups: dict[tuple[str, str], list[RebaseFamilyRescue]] = {}
    for rescue in rescues:
        key = (rescue.motif_iupac, rescue.methylation)
        groups.setdefault(key, []).append(rescue)

    ranked = sorted(
        groups.items(),
        key=lambda item: (
            -len(item[1]),
            -sum(rescue.bits for rescue in item[1]),
            item[0][0],
            item[0][1],
        ),
    )
    best_key, best_group = ranked[0]
    if len(best_group) < 2:
        return None

    representative = sorted(
        best_group,
        key=lambda rescue: (
            rescue.evalue,
            -rescue.bits,
            -rescue.pident,
            -rescue.qcov,
            -rescue.tcov,
        ),
    )[0]
    return RebaseFamilyRescue(
        source_hit=representative.source_hit,
        mapped_subject_id=representative.mapped_subject_id,
        motif_iupac=best_key[0],
        methylation=best_key[1],
        evalue=representative.evalue,
        bits=representative.bits,
        pident=representative.pident,
        qcov=representative.qcov,
        tcov=representative.tcov,
    )


def choose_consensus_rebase_hint(
    hints: list[RebaseHint],
) -> Optional[RebaseHint]:
    if not hints:
        return None

    groups: dict[tuple[str, str], list[RebaseHint]] = {}
    for hint in hints:
        key = (hint.motif_iupac, normalize_methylation(hint.methylation))
        groups.setdefault(key, []).append(hint)

    ranked = sorted(
        groups.items(),
        key=lambda item: (
            -len(item[1]),
            -sum(hint.bits for hint in item[1]),
            item[0][0],
            item[0][1],
        ),
    )
    best_key, best_group = ranked[0]
    representative = sorted(
        best_group,
        key=lambda hint: (
            hint.evalue,
            -hint.bits,
            -hint.pident,
            -hint.qcov,
            -hint.tcov,
        ),
    )[0]
    return RebaseHint(
        source_hit=representative.source_hit,
        mapped_subject_id=representative.mapped_subject_id,
        motif_iupac=best_key[0],
        methylation=best_key[1],
        evalue=representative.evalue,
        bits=representative.bits,
        pident=representative.pident,
        qcov=representative.qcov,
        tcov=representative.tcov,
    )


def confidence_label(hit: Hit) -> str:
    if hit.evalue <= 1e-50 and hit.pident >= 50 and hit.qcov >= 80 and hit.tcov >= 80:
        return "high"
    if hit.evalue <= 1e-20 and hit.pident >= 30 and hit.qcov >= 70 and hit.tcov >= 70:
        return "medium"
    return "low"


def confidence_from_rebase_family_rescue(rescue: RebaseFamilyRescue) -> str:
    if rescue.pident >= 65.0 and rescue.qcov >= 90.0 and rescue.tcov >= 90.0:
        return "high"
    if rescue.pident >= 50.0 and rescue.qcov >= 80.0 and rescue.tcov >= 80.0:
        return "medium"
    return "low"


def confidence_from_rebase_hint(hint: RebaseHint) -> str:
    if hint.bits >= 90.0 and hint.pident >= 28.0 and hint.qcov >= 90.0:
        return "medium"
    return "low"


def max_cluster_confidence_score(confidence: str) -> int:
    value = confidence.strip().lower()
    if value == "high":
        return 2
    if value == "medium":
        return 1
    return 0


def infer_named_rebase_motif(target_id: str) -> Optional[tuple[str, str]]:
    name = target_id.lower()
    if "dcm" in name:
        return ("CCWGG", "m5C")
    if "dam" in name:
        return ("GATC", "m6A")
    return None


def _percent_if_fraction(value: float) -> float:
    if 0.0 <= value <= 1.0:
        return value * 100.0
    return value


def _optional_float(value: str) -> Optional[float]:
    text = value.strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _format_optional_float(value: Optional[float]) -> str:
    if value is None:
        return ""
    return f"{value:g}"


def _safe_int(value: str) -> int:
    text = value.strip()
    if not text:
        return 0
    try:
        return int(float(text))
    except ValueError:
        return 0
