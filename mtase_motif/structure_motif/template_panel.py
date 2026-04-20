from __future__ import annotations

import shutil
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

from mtase_motif.mmcif import extract_dna_entity_sequences
from mtase_motif.motifs import normalize_iupac, pwm_to_iupac, sequences_to_pwm, write_meme_matrix
from mtase_motif.panel import add_weighted_kmers, choose_k_by_enrichment, extract_pdb_id, select_top_kmers
from mtase_motif.util import MtaseMotifError

from .foldseek import FoldseekHit


@dataclass(frozen=True)
class PanelBuild:
    motif_iupac: str
    pwm: list[list[float]]
    k: int
    nsites: int


@dataclass(frozen=True)
class TemplatePanelMotif:
    motif_iupac: str
    pwm: list[list[float]]
    k: int
    nsites: int
    best_target: str
    best_evalue: float
    best_bits: float
    best_alntmscore: float
    templates_used: int
    contributing_targets: tuple[str, ...]
    confidence: str


def build_panel_from_dna_sequences(
    dna_sequences: Sequence[str],
    *,
    weights: Optional[Sequence[float]] = None,
    kmin: int,
    kmax: int,
    panel_size: int,
) -> PanelBuild | None:
    clean_sequences = [seq.strip().upper() for seq in dna_sequences if seq and seq.strip()]
    if not clean_sequences:
        return None
    if weights is not None and len(weights) != len(clean_sequences):
        raise MtaseMotifError("weights length must match dna_sequences length")

    counts_by_k: dict[int, dict[str, float]] = {k: {} for k in range(int(kmin), int(kmax) + 1)}
    for idx, seq in enumerate(clean_sequences):
        weight = float(weights[idx]) if weights is not None else 1.0
        for k, counts in counts_by_k.items():
            add_weighted_kmers(counts, seq, k=k, weight=weight, canonicalize=True)

    best_k = choose_k_by_enrichment(counts_by_k)
    if best_k is None:
        return None
    best_counts = counts_by_k.get(best_k, {})
    kmers, kmer_weights = select_top_kmers(best_counts, limit=int(panel_size))
    if not kmers:
        return None

    pwm = sequences_to_pwm(kmers, weights=kmer_weights)
    motif_iupac = normalize_iupac(pwm_to_iupac(pwm))
    if not motif_iupac:
        return None

    return PanelBuild(
        motif_iupac=motif_iupac,
        pwm=pwm,
        k=int(best_k),
        nsites=len(kmers),
    )


def infer_template_panel_motif(
    hits: list[FoldseekHit],
    *,
    mmcif_mirror_dir: Optional[Path],
    mmcif_cache_dir: Path,
    allow_download: bool,
    top_hits: int,
    min_alntmscore: float,
    kmin: int,
    kmax: int,
    panel_size: int,
) -> TemplatePanelMotif | None:
    if not hits:
        return None
    mmcif_cache_dir.mkdir(parents=True, exist_ok=True)

    usable_hits = [h for h in hits if h.alntmscore >= min_alntmscore]
    if not usable_hits:
        return None
    usable_hits.sort(key=lambda h: (-h.alntmscore, h.evalue, -h.bits, h.target))
    selected = usable_hits[: max(1, top_hits)]

    weighted_sequences: list[str] = []
    weights: list[float] = []
    contributing_hits: list[FoldseekHit] = []
    for hit in selected:
        pdb_id = extract_pdb_id(hit.target)
        if pdb_id is None:
            continue
        cif = resolve_mmcif_path(
            pdb_id,
            mmcif_mirror_dir=mmcif_mirror_dir,
            mmcif_cache_dir=mmcif_cache_dir,
            allow_download=allow_download,
        )
        if cif is None:
            continue
        dna = extract_dna_entity_sequences(cif)
        if not dna:
            continue
        contributing_hits.append(hit)
        for seq in dna:
            weighted_sequences.append(seq)
            weights.append(hit.alntmscore)

    panel = build_panel_from_dna_sequences(
        weighted_sequences,
        weights=weights,
        kmin=kmin,
        kmax=kmax,
        panel_size=panel_size,
    )
    if panel is None:
        return None

    best = _best_contributing_hit(contributing_hits)
    if best is None:
        return None
    return TemplatePanelMotif(
        motif_iupac=panel.motif_iupac,
        pwm=panel.pwm,
        k=panel.k,
        nsites=panel.nsites,
        best_target=best.target,
        best_evalue=best.evalue,
        best_bits=best.bits,
        best_alntmscore=best.alntmscore,
        templates_used=len(contributing_hits),
        contributing_targets=tuple(hit.target for hit in contributing_hits),
        confidence=confidence_from_template_panel(best.alntmscore, len(contributing_hits)),
    )


def confidence_from_template_panel(best_alntmscore: float, templates_used: int) -> str:
    if best_alntmscore >= 0.7 and templates_used >= 3:
        return "high"
    if best_alntmscore >= 0.6 and templates_used >= 2:
        return "medium"
    return "low"


def _best_contributing_hit(hits: list[FoldseekHit]) -> FoldseekHit | None:
    if not hits:
        return None
    return sorted(hits, key=lambda h: (-h.alntmscore, h.evalue, -h.bits, h.target))[0]


def write_panel_outputs(
    consensus_path: Path,
    pwm_path: Path,
    candidate_id: str,
    panel: TemplatePanelMotif,
) -> None:
    consensus_path.parent.mkdir(parents=True, exist_ok=True)
    consensus_path.write_text(panel.motif_iupac + "\n")
    write_meme_matrix(pwm_path, motif_id=candidate_id, matrix=panel.pwm, nsites=panel.nsites)


def resolve_mmcif_path(
    pdb_id: str,
    *,
    mmcif_mirror_dir: Optional[Path],
    mmcif_cache_dir: Path,
    allow_download: bool,
) -> Optional[Path]:
    pid = pdb_id.lower()
    for root in [mmcif_cache_dir, mmcif_mirror_dir]:
        if root is None:
            continue
        found = find_existing_mmcif(root, pid)
        if found is not None:
            return found
    if not allow_download:
        return None

    dest = mmcif_cache_dir / f"{pid}.cif.gz"
    if dest.exists() and dest.stat().st_size > 0:
        return dest

    url = f"https://files.rcsb.org/download/{pid.upper()}.cif.gz"
    try:
        _download_url(url, dest)
        return dest
    except Exception:
        dest.unlink(missing_ok=True)
        return None


def find_existing_mmcif(root: Path, pdb_id: str) -> Optional[Path]:
    pid = pdb_id.lower()
    sub = pid[1:3]
    candidates = [
        root / f"{pid}.cif.gz",
        root / f"{pid}.cif",
        root / f"{pid}.mmcif.gz",
        root / f"{pid}.mmcif",
        root / sub / f"{pid}.cif.gz",
        root / sub / f"{pid}.cif",
        root / sub / f"{pid}.mmcif.gz",
        root / sub / f"{pid}.mmcif",
    ]
    for path in candidates:
        if path.exists() and path.stat().st_size > 0:
            return path
    return None


def _download_url(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "mtase-motif/0.1"})
    with urllib.request.urlopen(req, timeout=90) as response, dest.open("wb") as out:
        shutil.copyfileobj(response, out)
    if not dest.exists() or dest.stat().st_size == 0:
        raise MtaseMotifError(f"Downloaded empty file from {url}")
