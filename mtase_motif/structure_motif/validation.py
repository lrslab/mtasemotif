from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

from mtase_motif.methylation import methylation_labels, normalize_methylation
from mtase_motif.motifs import IUPAC_TO_BASES, normalize_iupac, reverse_complement_iupac


@dataclass(frozen=True)
class StructureValidationResult:
    candidate_id: str
    gold_motif_iupac: str
    predicted_motif_iupac: str
    gold_methylation: str
    predicted_methylation: str
    match_type: str
    orientation: str
    positional_overlap: float
    methylation_match: bool
    method: str
    confidence: str


@dataclass(frozen=True)
class StructureValidationSummary:
    total_gold: int
    resolved_predictions: int
    exact_matches: int
    reverse_complement_matches: int
    compatible_matches: int
    methylation_matches: int
    unresolved_predictions: int
    mean_positional_overlap: float


def compare_iupac_motifs(gold_motif: str, predicted_motif: str) -> tuple[str, str, float]:
    gold = normalize_iupac(gold_motif)
    predicted = normalize_iupac(predicted_motif)
    if not gold or not predicted:
        return "missing_prediction", "forward", 0.0

    forward = _compare_oriented_motifs(gold, predicted)
    reverse = _compare_oriented_motifs(gold, reverse_complement_iupac(predicted))
    if reverse > forward:
        orientation = "reverse"
        best = reverse
    else:
        orientation = "forward"
        best = forward

    exact, compatible, overlap, length_match = best
    if exact and orientation == "forward":
        return "exact", orientation, overlap
    if exact:
        return "reverse_complement", orientation, overlap
    if compatible:
        return "compatible", orientation, overlap
    if not length_match:
        return "length_mismatch", orientation, overlap
    return "incompatible", orientation, overlap


def validate_structure_predictions(
    gold_tsv: Path,
    predictions_tsv: Path,
    out_dir: Optional[Path] = None,
    *,
    split_filter: Optional[set[str]] = None,
) -> StructureValidationSummary:
    gold_rows = _read_tsv_rows(Path(gold_tsv))
    pred_rows = _read_tsv_rows(Path(predictions_tsv))
    predictions = {row.get("candidate_id", "").strip(): row for row in pred_rows if row.get("candidate_id", "").strip()}

    results: list[StructureValidationResult] = []
    for row in gold_rows:
        candidate_id = row.get("candidate_id", "").strip() or row.get("target", "").strip()
        if not candidate_id:
            continue
        split = row.get("split", "").strip()
        if split_filter and split not in split_filter:
            continue
        gold_motif = normalize_iupac(row.get("motif_iupac", ""))
        gold_methylation = normalize_methylation(row.get("methylation", ""))
        pred = predictions.get(candidate_id, {})
        predicted_motif = normalize_iupac(pred.get("motif_iupac", ""))
        predicted_methylation = normalize_methylation(pred.get("methylation", ""))
        match_type, orientation, overlap = compare_iupac_motifs(gold_motif, predicted_motif)
        results.append(
            StructureValidationResult(
                candidate_id=candidate_id,
                gold_motif_iupac=gold_motif,
                predicted_motif_iupac=predicted_motif,
                gold_methylation=gold_methylation,
                predicted_methylation=predicted_methylation,
                match_type=match_type,
                orientation=orientation,
                positional_overlap=overlap,
                methylation_match=_methylation_matches(gold_methylation, predicted_methylation),
                method=pred.get("method", "").strip(),
                confidence=pred.get("confidence", "").strip(),
            )
        )

    summary = _summarize_validation_results(results)
    if out_dir is not None:
        _write_validation_outputs(Path(out_dir), results, summary)
    return summary


def _compare_oriented_motifs(gold: str, predicted: str) -> tuple[bool, bool, float, bool]:
    if len(gold) != len(predicted):
        return False, False, 0.0, False

    exact = True
    compatible = True
    overlaps: list[float] = []
    for gold_ch, pred_ch in zip(gold, predicted):
        gold_bases = IUPAC_TO_BASES.get(gold_ch, IUPAC_TO_BASES["N"])
        pred_bases = IUPAC_TO_BASES.get(pred_ch, IUPAC_TO_BASES["N"])
        if gold_bases != pred_bases:
            exact = False
        intersection = gold_bases & pred_bases
        if not intersection:
            compatible = False
            overlaps.append(0.0)
            continue
        union = gold_bases | pred_bases
        overlaps.append(len(intersection) / float(len(union)))
    overlap = sum(overlaps) / float(len(overlaps)) if overlaps else 0.0
    return exact, compatible, overlap, True


def _methylation_matches(gold: str, predicted: str) -> bool:
    gold_norm = normalize_methylation(gold)
    predicted_norm = normalize_methylation(predicted)
    if not gold_norm or not predicted_norm:
        return False
    if gold_norm == predicted_norm:
        return True
    return bool(methylation_labels(gold_norm) & methylation_labels(predicted_norm))


def _summarize_validation_results(results: list[StructureValidationResult]) -> StructureValidationSummary:
    resolved = [result for result in results if result.predicted_motif_iupac]
    overlaps = [result.positional_overlap for result in resolved]
    return StructureValidationSummary(
        total_gold=len(results),
        resolved_predictions=len(resolved),
        exact_matches=sum(1 for result in results if result.match_type == "exact"),
        reverse_complement_matches=sum(1 for result in results if result.match_type == "reverse_complement"),
        compatible_matches=sum(1 for result in results if result.match_type in {"exact", "reverse_complement", "compatible"}),
        methylation_matches=sum(1 for result in results if result.methylation_match),
        unresolved_predictions=sum(1 for result in results if result.match_type == "missing_prediction"),
        mean_positional_overlap=(sum(overlaps) / float(len(overlaps))) if overlaps else 0.0,
    )


def _write_validation_outputs(
    out_dir: Path,
    results: list[StructureValidationResult],
    summary: StructureValidationSummary,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    results_path = out_dir / "validation_results.tsv"
    with results_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = list(asdict(results[0]).keys()) if results else list(asdict(_empty_result()).keys())
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for result in results:
            writer.writerow(asdict(result))
    (out_dir / "validation_summary.json").write_text(json.dumps(asdict(summary), indent=2, sort_keys=True))


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        return [{key: (value or "") for key, value in row.items()} for row in reader]


def _empty_result() -> StructureValidationResult:
    return StructureValidationResult(
        candidate_id="",
        gold_motif_iupac="",
        predicted_motif_iupac="",
        gold_methylation="",
        predicted_methylation="",
        match_type="missing_prediction",
        orientation="forward",
        positional_overlap=0.0,
        methylation_match=False,
        method="",
        confidence="",
    )
