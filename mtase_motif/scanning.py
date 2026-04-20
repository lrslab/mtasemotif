from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from mtase_motif.motifs import (
    canonicalize_iupac_orientation,
    classify_iupac_motif,
    infer_mod_position,
    is_palindromic_iupac,
    normalize_iupac,
    reverse_complement_iupac,
)
from mtase_motif.qc import can_use_exact_iupac_fallback, parse_fimo_hits, write_exact_iupac_fimo
from mtase_motif.util import MtaseMotifError


LOW_CANDIDATE_SUPPORT_THRESHOLD = 0.2
SUMMARY_COLS = [
    "candidate_id",
    "contig",
    "start",
    "end",
    "strand",
    "type_hint",
    "candidate_methylation_hint",
    "candidate_confidence",
    "predicted_motif",
    "predicted_motif_canonical",
    "predicted_motif_reverse_complement",
    "mod_position",
    "motif_class",
    "related_motif",
    "related_candidate_id",
    "related_confidence",
    "related_method",
    "hint_motif",
    "hint_methylation",
    "hint_hit_id",
    "hint_confidence",
    "hint_method",
    "assignment_state",
    "confidence",
    "num_sites",
    "density_per_mb",
    "warnings",
    "evidence_method",
]


def scan_and_summarize(
    *,
    genome_fa: Path,
    candidates_tsv: Path,
    out_summary: Path,
    candidate_id: str = "",
    call_tsv: Optional[Path] = None,
    motifs_tsv: Optional[Path] = None,
    pwm_path: Optional[Path] = None,
    out_fimo: Optional[Path] = None,
    out_qc: Optional[Path] = None,
    fimo_thresh: float = 1e-4,
    bg_order: int = 0,
) -> list[dict[str, str]]:
    candidate_mode = bool(candidate_id)
    if candidate_mode:
        if call_tsv is None or pwm_path is None or out_fimo is None or out_qc is None:
            raise MtaseMotifError(
                "candidate_id mode requires call_tsv, pwm_path, out_fimo, and out_qc"
            )
    else:
        if motifs_tsv is None:
            raise MtaseMotifError("motifs_tsv is required when candidate_id is not set")

    out_dir = out_summary.parent
    work_dir = out_dir / "work" / "fimo"
    work_dir.mkdir(parents=True, exist_ok=True)

    bg_file = work_dir / "bg.txt"
    genome_len = genome_length(genome_fa)
    if bg_order == 0:
        write_bgfile(bg_file, genome_base_freqs(genome_fa))
    else:
        fasta_get_markov = shutil.which("fasta-get-markov")
        if fasta_get_markov is None:
            raise MtaseMotifError(
                "bg_order>0 requested but fasta-get-markov not found in PATH (MEME suite)."
            )
        proc = subprocess.run(
            [fasta_get_markov, "-m", str(bg_order), str(genome_fa), str(bg_file)],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            raise MtaseMotifError(proc.stderr.strip() or "fasta-get-markov failed")

    candidates = load_candidates(candidates_tsv)
    motif_calls = load_motif_calls(call_tsv if candidate_mode else motifs_tsv)

    fimo = shutil.which("fimo")
    if fimo is None:
        raise MtaseMotifError("fimo not found in PATH (MEME suite required for scanning)")

    rows: list[dict[str, str]] = []
    candidate_ids = [candidate_id] if candidate_mode else list(candidates.keys())
    for cid in candidate_ids:
        candidate = candidates.get(cid)
        if candidate is None:
            continue
        motif = motif_calls.get(cid, {})
        motif_iupac = normalize_iupac(motif.get("motif_iupac", ""))
        hint_motif = normalize_iupac(motif.get("hint_motif_iupac", ""))
        evidence_method = motif.get("method", "")
        motif_conf = motif.get("confidence", "")
        type_hint = candidate.get("type_hint", "")
        candidate_methylation_hint = candidate.get("methylation_hint", "")
        candidate_confidence = candidate.get("confidence", "")
        predicted_motif_canonical = motif.get("motif_canonical_iupac", "") or canonicalize_iupac_orientation(
            motif_iupac
        )
        predicted_motif_rc = motif.get("motif_reverse_complement_iupac", "") or (
            reverse_complement_iupac(motif_iupac) if motif_iupac else ""
        )
        mod_position = motif.get("mod_position", "")
        if not mod_position:
            inferred_mod_position = infer_mod_position(motif_iupac, motif.get("methylation", ""))
            mod_position = str(inferred_mod_position) if inferred_mod_position is not None else ""
        motif_class = motif.get("motif_class", "") or classify_iupac_motif(motif_iupac)
        assignment_state = motif.get("assignment_state", "") or infer_assignment_state(motif)

        warnings: list[str] = []
        num_sites = 0
        density = 0.0
        strand_balance = 0.0
        is_palindromic = False

        if motif_iupac:
            candidate_pwm = pwm_path if candidate_mode else out_dir / cid / "motif" / "pwm.meme"
            candidate_fimo = out_fimo if candidate_mode else out_dir / cid / "fimo" / "fimo.tsv"
            run_fimo(fimo, candidate_pwm, genome_fa, bg_file, candidate_fimo.parent, fimo_thresh)
            num_sites, density, strand_balance = parse_fimo_hits(candidate_fimo, genome_len)
            if num_sites == 0 and can_use_exact_iupac_fallback(motif_iupac):
                write_exact_iupac_fimo(
                    candidate_fimo,
                    genome_fa=genome_fa,
                    motif_id=cid,
                    motif_iupac=motif_iupac,
                )
                num_sites, density, strand_balance = parse_fimo_hits(candidate_fimo, genome_len)

            if num_sites == 0:
                warnings.append("no_sites")
            if density > 5000:
                warnings.append("too_many_sites")
            if 0 < num_sites and strand_balance > 0.8:
                warnings.append("strand_imbalance")

            is_palindromic = is_palindromic_iupac(motif_iupac)
            if is_type_ii_like(type_hint) and not is_palindromic:
                warnings.append("type_ii_non_palindromic")

            qc_path = out_qc if candidate_mode else out_dir / cid / "qc" / "qc.json"
            write_qc_json(
                qc_path,
                num_sites=num_sites,
                density=density,
                strand_balance=strand_balance,
                motif_iupac=motif_iupac,
                is_palindromic=is_palindromic,
                warnings=warnings,
            )
        else:
            warnings.append("no_motif")
            if hint_motif:
                warnings.append("hint_only")
            if is_low_support_candidate(candidate_confidence):
                warnings.append("low_candidate_support")
            candidate_fimo = out_fimo if candidate_mode else out_dir / cid / "fimo" / "fimo.tsv"
            write_empty_fimo(candidate_fimo)
            qc_path = out_qc if candidate_mode else out_dir / cid / "qc" / "qc.json"
            write_qc_json(
                qc_path,
                num_sites=0,
                density=0.0,
                strand_balance=0.0,
                motif_iupac="",
                is_palindromic=False,
                warnings=warnings,
            )

        rows.append(
            {
                "candidate_id": cid,
                "contig": candidate.get("contig", ""),
                "start": candidate.get("start", ""),
                "end": candidate.get("end", ""),
                "strand": candidate.get("strand", ""),
                "type_hint": type_hint,
                "candidate_methylation_hint": candidate_methylation_hint,
                "candidate_confidence": candidate_confidence,
                "predicted_motif": motif_iupac,
                "predicted_motif_canonical": predicted_motif_canonical,
                "predicted_motif_reverse_complement": predicted_motif_rc,
                "mod_position": mod_position,
                "motif_class": motif_class,
                "related_motif": motif.get("related_motif_iupac", ""),
                "related_candidate_id": motif.get("related_candidate_id", ""),
                "related_confidence": motif.get("related_confidence", ""),
                "related_method": motif.get("related_method", ""),
                "hint_motif": hint_motif,
                "hint_methylation": motif.get("hint_methylation", ""),
                "hint_hit_id": motif.get("hint_hit_id", ""),
                "hint_confidence": motif.get("hint_confidence", ""),
                "hint_method": motif.get("hint_method", ""),
                "assignment_state": assignment_state,
                "confidence": motif_conf,
                "num_sites": str(num_sites),
                "density_per_mb": f"{density:.3f}",
                "warnings": ";".join(warnings),
                "evidence_method": evidence_method,
            }
        )

    rows.sort(key=lambda row: (row.get("contig") or "~", _safe_int(row.get("start")), row["candidate_id"]))
    write_summary_rows(out_summary, SUMMARY_COLS, rows)
    return rows


def genome_length(genome_fa: Path) -> int:
    total_len = 0
    with genome_fa.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith(">"):
                continue
            total_len += len(line)
    return total_len


def genome_base_freqs(genome_fa: Path) -> dict[str, float]:
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    with genome_fa.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith(">"):
                continue
            for base in line.upper():
                if base in counts:
                    counts[base] += 1
    total_bases = sum(counts.values())
    if total_bases == 0:
        return {base: 0.25 for base in "ACGT"}
    return {base: counts[base] / total_bases for base in "ACGT"}


def write_bgfile(path: Path, freqs: dict[str, float]) -> None:
    path.write_text(
        "# 0-order Markov background model\n"
        + "\n".join([f"{base} {freqs[base]:.6f}" for base in "ACGT"])
        + "\n"
    )


def load_candidates(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    lines = path.read_text().splitlines()
    if not lines:
        return rows
    header = lines[0].split("\t")
    if "candidate_id" not in header:
        raise MtaseMotifError(f"Malformed candidates TSV (missing candidate_id column): {path}")
    for raw in lines[1:]:
        if not raw:
            continue
        parts = raw.split("\t")
        row = {header[i]: parts[i] if i < len(parts) else "" for i in range(len(header))}
        candidate_id = row.get("candidate_id")
        if candidate_id:
            rows[candidate_id] = row
    return rows


def load_motif_calls(path: Optional[Path]) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    if path is None:
        return rows
    lines = path.read_text().splitlines()
    if not lines:
        return rows
    header = lines[0].split("\t")
    if "candidate_id" not in header:
        raise MtaseMotifError(f"Malformed motif calls TSV (missing candidate_id column): {path}")
    for raw in lines[1:]:
        if not raw:
            continue
        parts = raw.split("\t")
        row = {header[i]: parts[i] if i < len(parts) else "" for i in range(len(header))}
        candidate_id = row.get("candidate_id")
        if candidate_id:
            rows[candidate_id] = row
    return rows


def infer_assignment_state(motif_row: dict[str, str]) -> str:
    if normalize_iupac(motif_row.get("motif_iupac", "")):
        return "linked"
    if normalize_iupac(motif_row.get("related_motif_iupac", "")) or normalize_iupac(
        motif_row.get("hint_motif_iupac", "")
    ):
        return "ambiguous"
    return "unresolved"


def is_low_support_candidate(confidence_raw: str) -> bool:
    try:
        confidence = float(confidence_raw)
    except Exception:
        return False
    return confidence < LOW_CANDIDATE_SUPPORT_THRESHOLD


def run_fimo(
    fimo_exe: str,
    pwm_meme: Path,
    genome_fa: Path,
    bg_file: Path,
    out_dir: Path,
    thresh: float,
) -> None:
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        [
            fimo_exe,
            "--verbosity",
            "1",
            "--thresh",
            str(thresh),
            "--bgfile",
            str(bg_file),
            "--oc",
            str(out_dir),
            str(pwm_meme),
            str(genome_fa),
        ],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise MtaseMotifError(proc.stderr.strip() or "fimo failed")


def write_empty_fimo(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n"
    )


def write_qc_json(
    path: Path,
    *,
    num_sites: int,
    density: float,
    strand_balance: float,
    motif_iupac: str,
    is_palindromic: bool,
    warnings: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "num_sites": num_sites,
                "density_per_mb": density,
                "strand_balance": strand_balance,
                "motif_iupac": motif_iupac,
                "is_palindromic": is_palindromic,
                "warnings": warnings,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def write_summary_rows(path: Path, cols: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("\t".join(cols) + "\n")
        for row in rows:
            handle.write("\t".join(row.get(col, "") for col in cols) + "\n")


def _safe_int(value: Optional[str]) -> int:
    try:
        return int(value or "0")
    except Exception:
        return 0


def is_type_ii_like(type_hint: str) -> bool:
    value = (type_hint or "").lower()
    return value.startswith("type_ii") or value.startswith("orphan_like")
