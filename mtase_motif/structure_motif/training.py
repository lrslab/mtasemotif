from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

from mtase_motif.methylation import normalize_methylation
from mtase_motif.mmcif import extract_dna_entity_sequences
from mtase_motif.motifs import iupac_degeneracy, normalize_iupac
from mtase_motif.panel import extract_pdb_id

from .template_panel import build_panel_from_dna_sequences, resolve_mmcif_path


@dataclass(frozen=True)
class StructureTrainingRecord:
    target: str
    motif_iupac: str
    methylation: str
    family_id: str
    split: str
    pdb_id: str
    structure_path: str
    has_structure: bool
    has_dna: bool
    dna_entities: int
    usable_for_label_transfer: bool
    usable_for_template_panel: bool
    self_panel_iupac: str
    self_panel_k: int
    motif_width: int
    motif_degeneracy: int


@dataclass(frozen=True)
class StructureTrainingSummary:
    total_rows: int
    label_transfer_ready: int
    template_panel_ready: int
    missing_structure: int
    missing_dna: int
    unique_motifs: int
    unique_families: int
    split_counts: dict[str, int]


def build_structure_label_db(
    manifest_tsv: Path,
    out_dir: Path,
    *,
    mmcif_mirror_dir: Optional[Path] = None,
    mmcif_cache_dir: Optional[Path] = None,
    allow_download: bool = False,
    kmin: int = 6,
    kmax: int = 8,
    panel_size: int = 25,
) -> StructureTrainingSummary:
    manifest_tsv = Path(manifest_tsv)
    out_dir = Path(out_dir)
    mmcif_cache = mmcif_cache_dir or (out_dir / "mmcif_cache")
    rows = _read_tsv_rows(manifest_tsv)
    records = [
        _build_training_record(
            row,
            manifest_dir=manifest_tsv.parent,
            mmcif_mirror_dir=mmcif_mirror_dir,
            mmcif_cache_dir=mmcif_cache,
            allow_download=allow_download,
            kmin=kmin,
            kmax=kmax,
            panel_size=panel_size,
        )
        for row in rows
    ]

    out_dir.mkdir(parents=True, exist_ok=True)
    _write_labels_tsv(out_dir / "labels.tsv", records)
    _write_training_examples_tsv(out_dir / "training_examples.tsv", records)
    summary = _summarize_training_records(records)
    (out_dir / "training_summary.json").write_text(json.dumps(asdict(summary), indent=2, sort_keys=True))
    return summary


def _build_training_record(
    row: dict[str, str],
    *,
    manifest_dir: Path,
    mmcif_mirror_dir: Optional[Path],
    mmcif_cache_dir: Path,
    allow_download: bool,
    kmin: int,
    kmax: int,
    panel_size: int,
) -> StructureTrainingRecord:
    target = row.get("target", "").strip()
    motif_iupac = normalize_iupac(row.get("motif_iupac", ""))
    methylation = normalize_methylation(row.get("methylation", ""))
    family_id = row.get("family_id", "").strip()
    split = row.get("split", "").strip()
    pdb_id = row.get("pdb_id", "").strip().lower() or (extract_pdb_id(target) or "")

    structure_path = _resolve_structure_path(
        row,
        manifest_dir=manifest_dir,
        mmcif_mirror_dir=mmcif_mirror_dir,
        mmcif_cache_dir=mmcif_cache_dir,
        allow_download=allow_download,
        pdb_id=pdb_id,
    )
    dna_sequences: list[str] = []
    if structure_path is not None and structure_path.suffix in {".cif", ".mmcif", ".gz"}:
        try:
            dna_sequences = extract_dna_entity_sequences(structure_path)
        except Exception:
            dna_sequences = []

    panel = build_panel_from_dna_sequences(
        dna_sequences,
        kmin=kmin,
        kmax=kmax,
        panel_size=panel_size,
    )
    return StructureTrainingRecord(
        target=target,
        motif_iupac=motif_iupac,
        methylation=methylation,
        family_id=family_id,
        split=split,
        pdb_id=pdb_id,
        structure_path=str(structure_path) if structure_path is not None else "",
        has_structure=structure_path is not None,
        has_dna=bool(dna_sequences),
        dna_entities=len(dna_sequences),
        usable_for_label_transfer=bool(target and motif_iupac),
        usable_for_template_panel=bool(motif_iupac and panel is not None),
        self_panel_iupac=panel.motif_iupac if panel is not None else "",
        self_panel_k=panel.k if panel is not None else 0,
        motif_width=len(motif_iupac),
        motif_degeneracy=iupac_degeneracy(motif_iupac) if motif_iupac else 0,
    )


def _resolve_structure_path(
    row: dict[str, str],
    *,
    manifest_dir: Path,
    mmcif_mirror_dir: Optional[Path],
    mmcif_cache_dir: Path,
    allow_download: bool,
    pdb_id: str,
) -> Optional[Path]:
    raw_path = row.get("structure_path", "").strip()
    if raw_path:
        candidate = Path(raw_path)
        if not candidate.is_absolute():
            candidate = manifest_dir / candidate
        if candidate.exists() and candidate.stat().st_size > 0:
            return candidate.resolve()
    if pdb_id:
        return resolve_mmcif_path(
            pdb_id,
            mmcif_mirror_dir=mmcif_mirror_dir,
            mmcif_cache_dir=mmcif_cache_dir,
            allow_download=allow_download,
        )
    return None


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        return [{key: (value or "") for key, value in row.items()} for row in reader]


def _write_labels_tsv(path: Path, records: list[StructureTrainingRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["target", "motif_iupac", "methylation"])
        for record in records:
            if not record.target or not record.motif_iupac:
                continue
            writer.writerow([record.target, record.motif_iupac, record.methylation])


def _write_training_examples_tsv(path: Path, records: list[StructureTrainingRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(asdict(records[0]).keys()) if records else list(asdict(_empty_record()).keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for record in records:
            writer.writerow(asdict(record))


def _summarize_training_records(records: list[StructureTrainingRecord]) -> StructureTrainingSummary:
    split_counts: dict[str, int] = {}
    motifs = {record.motif_iupac for record in records if record.motif_iupac}
    families = {record.family_id for record in records if record.family_id}
    for record in records:
        if record.split:
            split_counts[record.split] = split_counts.get(record.split, 0) + 1

    return StructureTrainingSummary(
        total_rows=len(records),
        label_transfer_ready=sum(1 for record in records if record.usable_for_label_transfer),
        template_panel_ready=sum(1 for record in records if record.usable_for_template_panel),
        missing_structure=sum(1 for record in records if not record.has_structure),
        missing_dna=sum(1 for record in records if record.has_structure and not record.has_dna),
        unique_motifs=len(motifs),
        unique_families=len(families),
        split_counts=split_counts,
    )


def _empty_record() -> StructureTrainingRecord:
    return StructureTrainingRecord(
        target="",
        motif_iupac="",
        methylation="",
        family_id="",
        split="",
        pdb_id="",
        structure_path="",
        has_structure=False,
        has_dna=False,
        dna_entities=0,
        usable_for_label_transfer=False,
        usable_for_template_panel=False,
        self_panel_iupac="",
        self_panel_k=0,
        motif_width=0,
        motif_degeneracy=0,
    )
