from __future__ import annotations

import logging
import math
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from mtase_motif.hmmer import DomainHit, parse_domtblout
from mtase_motif.util import MtaseMotifError


MTASE_CATALYTIC_PFAMS = {"PF01555", "PF00145"}
NEIGHBOR_WINDOW_BP = 10_000
REBASE_DAM_RESCUE_IDS = {
    "M.EcoKDam",
    "M.StyLT2Dam",
    "M.AacDam",
    "M.YpeDam",
    "M.CneDam",
}

CANDIDATE_COLS = [
    "candidate_id",
    "contig",
    "start",
    "end",
    "strand",
    "domains",
    "type_hint",
    "methylation_hint",
    "neighborhood_features",
    "confidence",
    "rationale",
]


log = logging.getLogger(__name__)


@dataclass(frozen=True)
class GeneCoord:
    contig: str
    start: int
    end: int
    strand: str


@dataclass(frozen=True)
class RoleFlags:
    has_mtase: bool
    has_type_i_m: bool
    has_type_i_s: bool
    has_type_i_r: bool
    has_type_iii_mod: bool
    has_type_iii_res: bool
    has_dam_like: bool
    has_dcm_like: bool
    has_likely_rease: bool


@dataclass(frozen=True)
class NeighborhoodSummary:
    type_i_s_ids: tuple[str, ...]
    type_i_r_ids: tuple[str, ...]
    type_iii_res_ids: tuple[str, ...]
    rease_ids: tuple[str, ...]
    min_type_i_s_dist: Optional[int]
    min_type_i_r_dist: Optional[int]
    min_type_iii_res_dist: Optional[int]
    min_rease_dist: Optional[int]

    @property
    def has_type_i_s(self) -> bool:
        return len(self.type_i_s_ids) > 0

    @property
    def has_type_i_r(self) -> bool:
        return len(self.type_i_r_ids) > 0

    @property
    def has_type_iii_res(self) -> bool:
        return len(self.type_iii_res_ids) > 0

    @property
    def has_likely_rease(self) -> bool:
        return len(self.rease_ids) > 0


def empty_neighborhood_summary() -> NeighborhoodSummary:
    return NeighborhoodSummary(
        type_i_s_ids=(),
        type_i_r_ids=(),
        type_iii_res_ids=(),
        rease_ids=(),
        min_type_i_s_dist=None,
        min_type_i_r_dist=None,
        min_type_iii_res_dist=None,
        min_rease_dist=None,
    )


@dataclass(frozen=True)
class RescueHit:
    target: str
    evalue: float
    bits: float
    pident: float
    qcov: float
    tcov: float


def build_candidates(
    *,
    pfam_domtbl_path: Path,
    tigr_domtbl_path: Path,
    gff_path: Path,
    out_path: Path,
    proteins_faa: Optional[Path] = None,
    rebase_proteins: Optional[Path] = None,
    pfam_evalue_thresh: float = 1e-5,
    tigr_evalue_thresh: float = 1e-5,
    max_domains: int = 5,
    rescue_max_evalue: float = 1e-20,
    rescue_min_pident: float = 70.0,
    rescue_min_qcov: float = 80.0,
    rescue_min_tcov: float = 80.0,
) -> list[dict[str, str]]:
    coords = parse_gff_coords(gff_path)

    hits_by_query: dict[str, list[DomainHit]] = {}
    _merge_hits(
        hits_by_query,
        parse_domtblout(pfam_domtbl_path, pfam_evalue_thresh, source="pfam", program="hmmscan"),
    )
    _merge_hits(
        hits_by_query,
        parse_domtblout(
            tigr_domtbl_path,
            tigr_evalue_thresh,
            source="tigrfams",
            program="hmmsearch",
        ),
    )

    rows: list[dict[str, str]] = []
    for query_id, hits in hits_by_query.items():
        hits_sorted = sorted(hits, key=lambda h: h.i_evalue)
        best = choose_best_mtase_hit(hits_sorted)
        if best is None:
            continue

        flags = infer_role_flags(hits_sorted)
        candidate_coord = coords.get(query_id)
        if candidate_coord is None:
            neighborhood = empty_neighborhood_summary()
        else:
            neighborhood = summarize_neighborhood(
                candidate_id=query_id,
                candidate_coord=candidate_coord,
                coords=coords,
                hits_by_query=hits_by_query,
                window_bp=NEIGHBOR_WINDOW_BP,
            )

        rows.append(
            {
                "candidate_id": query_id,
                **coord_fields(candidate_coord),
                "domains": ";".join(
                    [
                        f"{h.source}:{h.model_acc}({h.model_name}):{h.i_evalue:g}"
                        for h in hits_sorted[:max_domains]
                    ]
                ),
                "type_hint": classify_type_hint(best=best, flags=flags, neighborhood=neighborhood),
                "methylation_hint": methylation_hint_for(best=best, flags=flags),
                "neighborhood_features": neighborhood_features_text(
                    flags,
                    neighborhood,
                    coordinates_available=candidate_coord is not None,
                ),
                "confidence": f"{confidence_from_evalue(best.i_evalue):.3f}",
                "rationale": (
                    f"{best.source}_hit {best.model_acc}({best.model_name}) "
                    f"i_evalue={best.i_evalue:g}"
                ),
            }
        )

    rescue_hits = find_named_rebase_rescue_hits(
        proteins_faa=proteins_faa,
        rebase_proteins=rebase_proteins,
        max_evalue=rescue_max_evalue,
        min_pident=rescue_min_pident,
        min_qcov=rescue_min_qcov,
        min_tcov=rescue_min_tcov,
        work_dir=out_path.parent / "work" / "candidate_rescue",
    )
    existing_ids = {row["candidate_id"] for row in rows}
    for query_id, hit in rescue_hits.items():
        if query_id in existing_ids:
            continue
        candidate_coord = coords.get(query_id)
        methylation_hint = methylation_hint_for_named_rebase_target(hit.target)
        if not methylation_hint:
            continue
        rows.append(
            {
                "candidate_id": query_id,
                **coord_fields(candidate_coord),
                "domains": f"rebase_named:{hit.target}:{hit.evalue:g}",
                "type_hint": "type_II_or_orphan_like",
                "methylation_hint": methylation_hint,
                "neighborhood_features": (
                    "window_bp=10000;rebase_named_rescue=1"
                    if candidate_coord is not None
                    else "window_bp=10000;coord_source=unavailable;rebase_named_rescue=1"
                ),
                "confidence": f"{confidence_from_named_rebase_rescue(hit):.3f}",
                "rationale": (
                    f"rebase_named_rescue {hit.target} "
                    f"evalue={hit.evalue:g} pident={hit.pident:g} "
                    f"qcov={hit.qcov:g} tcov={hit.tcov:g}"
                ),
            }
        )

    rows.sort(key=_row_sort_key)
    write_candidates(out_path, rows)
    return rows


def write_candidates(out_path: Path, rows: list[dict[str, str]]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as handle:
        handle.write("\t".join(CANDIDATE_COLS) + "\n")
        for row in rows:
            handle.write("\t".join(row.get(col, "") for col in CANDIDATE_COLS) + "\n")


def _merge_hits(
    hits_by_query: dict[str, list[DomainHit]],
    new_hits: dict[str, list[DomainHit]],
) -> None:
    for query_id, hits in new_hits.items():
        hits_by_query.setdefault(query_id, []).extend(hits)


def choose_best_mtase_hit(hits_sorted: list[DomainHit]) -> Optional[DomainHit]:
    catalytic = [h for h in hits_sorted if h.source == "pfam" and h.model_acc in MTASE_CATALYTIC_PFAMS]
    if catalytic:
        return catalytic[0]

    tigr_mtase = [h for h in hits_sorted if h.source == "tigrfams" and _is_tigr_mtase_hit(h)]
    if tigr_mtase:
        return tigr_mtase[0]

    other_mtase = [h for h in hits_sorted if _is_methyltransferase_hit(h)]
    if other_mtase:
        return other_mtase[0]
    return None


def coord_fields(coord: Optional[GeneCoord]) -> dict[str, str]:
    if coord is None:
        return {
            "contig": "",
            "start": "",
            "end": "",
            "strand": "",
        }
    return {
        "contig": coord.contig,
        "start": str(coord.start),
        "end": str(coord.end),
        "strand": coord.strand,
    }


def find_named_rebase_rescue_hits(
    *,
    proteins_faa: Optional[Path],
    rebase_proteins: Optional[Path],
    max_evalue: float,
    min_pident: float,
    min_qcov: float,
    min_tcov: float,
    work_dir: Path,
) -> dict[str, RescueHit]:
    if proteins_faa is None or rebase_proteins is None:
        return {}
    if not proteins_faa.exists() or not rebase_proteins.exists():
        return {}

    blastp = shutil.which("blastp")
    if blastp is None:
        log.warning(
            "blastp not found in PATH; skipping named REBASE rescue against %s. "
            "Install BLAST+ to enable Dam/Dcm rescue candidates.",
            rebase_proteins,
        )
        return {}

    work_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="named_rebase_", dir=str(work_dir)) as tmpdir:
        subject_faa = Path(tmpdir) / "named_rebase_subset.faa"
        if write_named_rebase_subset(rebase_proteins, subject_faa) == 0:
            return {}

        proc = subprocess.run(
            [
                blastp,
                "-query",
                str(proteins_faa),
                "-subject",
                str(subject_faa),
                "-max_target_seqs",
                "10",
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
                    "blastp named REBASE rescue failed while searching "
                    f"{proteins_faa} against {subject_faa}"
                )
            )

    best_hits: dict[str, RescueHit] = {}
    for raw in proc.stdout.splitlines():
        if not raw.strip():
            continue
        parts = raw.split("\t")
        if len(parts) != 8:
            continue
        query_id, target_id, evalue_raw, bits_raw, pident_raw, qcov_raw, slen_raw, length_raw = parts
        try:
            evalue = float(evalue_raw)
            bits = float(bits_raw)
            pident = float(pident_raw)
            qcov = float(qcov_raw)
            slen = float(slen_raw)
            alen = float(length_raw)
        except ValueError:
            continue
        tcov = 0.0 if slen <= 0 else (100.0 * alen / slen)
        if evalue > max_evalue or pident < min_pident or qcov < min_qcov or tcov < min_tcov:
            continue
        hit = RescueHit(
            target=target_id,
            evalue=evalue,
            bits=bits,
            pident=pident,
            qcov=qcov,
            tcov=tcov,
        )
        current = best_hits.get(query_id)
        if current is None or (hit.evalue, -hit.bits, -hit.pident, -hit.qcov, -hit.tcov) < (
            current.evalue,
            -current.bits,
            -current.pident,
            -current.qcov,
            -current.tcov,
        ):
            best_hits[query_id] = hit
    return best_hits


def write_named_rebase_subset(rebase_proteins: Path, out_faa: Path) -> int:
    count = 0
    keep = False
    with rebase_proteins.open() as src, out_faa.open("w") as dst:
        for raw in src:
            if raw.startswith(">"):
                target_id = raw[1:].strip().split()[0]
                keep = target_id in REBASE_DAM_RESCUE_IDS
                if keep:
                    count += 1
                    dst.write(raw)
                continue
            if keep:
                dst.write(raw)
    return count


def methylation_hint_for_named_rebase_target(target_id: str) -> str:
    name = target_id.lower()
    if "dam" in name:
        return "m6A"
    if "dcm" in name:
        return "m5C"
    return ""


def confidence_from_named_rebase_rescue(hit: RescueHit) -> float:
    if hit.evalue <= 1e-50 and hit.pident >= 90 and hit.qcov >= 90 and hit.tcov >= 90:
        return 1.0
    if hit.evalue <= 1e-30 and hit.pident >= 80 and hit.qcov >= 85 and hit.tcov >= 85:
        return 0.9
    return 0.8


def infer_role_flags(hits: list[DomainHit]) -> RoleFlags:
    text_parts = [_hit_text(h) for h in hits]
    all_text = " || ".join(text_parts)
    has_mtase = any(_is_methyltransferase_hit(h) for h in hits)

    has_type_i_m = _has_any(
        all_text,
        ("hsdm", "type i methyltransferase", "type i m subunit", "m subunit"),
    ) and "type iii" not in all_text
    has_type_i_s = _has_any(
        all_text,
        ("hsds", "type i specificity", "specificity subunit", "type i s subunit"),
    )
    has_type_i_r = _has_any(
        all_text,
        ("hsdr", "type i restriction", "restriction subunit"),
    )
    has_type_iii_mod = _has_any(
        all_text,
        ("type iii methyltransferase", "type iii mod", "mod subunit"),
    )
    has_type_iii_res = _has_any(
        all_text,
        ("type iii restriction", "type iii res", "res subunit"),
    )
    has_dam_like = _has_any(
        all_text,
        (" dam", "dam ", "dna adenine methylase", "adenine methyltransferase"),
    )
    has_dcm_like = _has_any(
        all_text,
        (" dcm", "dcm ", "cytosine methyltransferase", "cytosine-specific"),
    )
    has_likely_rease = any(_is_likely_rease_hit(h) for h in hits)

    return RoleFlags(
        has_mtase=has_mtase,
        has_type_i_m=has_type_i_m,
        has_type_i_s=has_type_i_s,
        has_type_i_r=has_type_i_r,
        has_type_iii_mod=has_type_iii_mod,
        has_type_iii_res=has_type_iii_res,
        has_dam_like=has_dam_like,
        has_dcm_like=has_dcm_like,
        has_likely_rease=has_likely_rease,
    )


def summarize_neighborhood(
    *,
    candidate_id: str,
    candidate_coord: GeneCoord,
    coords: dict[str, GeneCoord],
    hits_by_query: dict[str, list[DomainHit]],
    window_bp: int,
) -> NeighborhoodSummary:
    type_i_s_dists: list[tuple[int, str]] = []
    type_i_r_dists: list[tuple[int, str]] = []
    type_iii_res_dists: list[tuple[int, str]] = []
    rease_dists: list[tuple[int, str]] = []

    for gene_id, gene_hits in hits_by_query.items():
        if gene_id == candidate_id:
            continue
        gene_coord = coords.get(gene_id)
        if gene_coord is None or gene_coord.contig != candidate_coord.contig:
            continue
        distance = gene_distance(candidate_coord, gene_coord)
        if distance > window_bp:
            continue

        gene_flags = infer_role_flags(gene_hits)
        if gene_flags.has_type_i_s:
            type_i_s_dists.append((distance, gene_id))
        if gene_flags.has_type_i_r:
            type_i_r_dists.append((distance, gene_id))
        if gene_flags.has_type_iii_res:
            type_iii_res_dists.append((distance, gene_id))
        if gene_flags.has_likely_rease:
            rease_dists.append((distance, gene_id))

    type_i_s_dists.sort()
    type_i_r_dists.sort()
    type_iii_res_dists.sort()
    rease_dists.sort()

    return NeighborhoodSummary(
        type_i_s_ids=tuple([gene_id for _, gene_id in type_i_s_dists]),
        type_i_r_ids=tuple([gene_id for _, gene_id in type_i_r_dists]),
        type_iii_res_ids=tuple([gene_id for _, gene_id in type_iii_res_dists]),
        rease_ids=tuple([gene_id for _, gene_id in rease_dists]),
        min_type_i_s_dist=type_i_s_dists[0][0] if type_i_s_dists else None,
        min_type_i_r_dist=type_i_r_dists[0][0] if type_i_r_dists else None,
        min_type_iii_res_dist=type_iii_res_dists[0][0] if type_iii_res_dists else None,
        min_rease_dist=rease_dists[0][0] if rease_dists else None,
    )


def gene_distance(a: GeneCoord, b: GeneCoord) -> int:
    if a.end < b.start:
        return b.start - a.end
    if b.end < a.start:
        return a.start - b.end
    return 0


def classify_type_hint(*, best: DomainHit, flags: RoleFlags, neighborhood: NeighborhoodSummary) -> str:
    if (flags.has_type_i_m or "hsdm" in _hit_text(best)) and neighborhood.has_type_i_s and neighborhood.has_type_i_r:
        return "type_I_like"
    if flags.has_type_iii_mod and neighborhood.has_type_iii_res:
        return "type_III_like"
    if (flags.has_dam_like or flags.has_dcm_like) and not neighborhood.has_likely_rease:
        return "orphan_like_dam_dcm"
    if neighborhood.has_likely_rease:
        return "type_II_rm_like"
    if best.source == "pfam" and best.model_acc in MTASE_CATALYTIC_PFAMS:
        return "type_II_or_orphan_like"
    if flags.has_mtase:
        return "mtase_like"
    return "unknown"


def neighborhood_features_text(
    flags: RoleFlags,
    neighborhood: NeighborhoodSummary,
    *,
    coordinates_available: bool = True,
) -> str:
    features: list[str] = [f"window_bp={NEIGHBOR_WINDOW_BP}"]
    if not coordinates_available:
        features.append("coord_source=unavailable")
    if flags.has_type_i_m:
        features.append("cand_typeI_M")
    if flags.has_type_iii_mod:
        features.append("cand_typeIII_Mod")
    if flags.has_dam_like:
        features.append("cand_dam_like")
    if flags.has_dcm_like:
        features.append("cand_dcm_like")
    if flags.has_likely_rease:
        features.append("cand_rease_like_domain")

    if neighborhood.type_i_s_ids:
        features.append("near_typeI_S=" + ",".join(neighborhood.type_i_s_ids[:3]))
    if neighborhood.type_i_r_ids:
        features.append("near_typeI_R=" + ",".join(neighborhood.type_i_r_ids[:3]))
    if neighborhood.type_iii_res_ids:
        features.append("near_typeIII_Res=" + ",".join(neighborhood.type_iii_res_ids[:3]))
    if neighborhood.rease_ids:
        features.append("near_rease=" + ",".join(neighborhood.rease_ids[:3]))
    if neighborhood.min_type_i_s_dist is not None:
        features.append(f"dist_typeI_S={neighborhood.min_type_i_s_dist}")
    if neighborhood.min_type_i_r_dist is not None:
        features.append(f"dist_typeI_R={neighborhood.min_type_i_r_dist}")
    if neighborhood.min_type_iii_res_dist is not None:
        features.append(f"dist_typeIII_Res={neighborhood.min_type_iii_res_dist}")
    if neighborhood.min_rease_dist is not None:
        features.append(f"dist_rease={neighborhood.min_rease_dist}")

    return ";".join(features)


def methylation_hint_for(*, best: DomainHit, flags: RoleFlags) -> str:
    if best.source == "pfam":
        if best.model_acc == "PF00145":
            return "m5C"
        if best.model_acc == "PF01555":
            return "m6A/m4C"
        return "unknown"

    if flags.has_dcm_like:
        return "m5C"
    if flags.has_dam_like:
        return "m6A"
    if flags.has_type_iii_mod:
        return "m6A/m4C"
    return "unknown"


def parse_gff_coords(gff: Path) -> dict[str, GeneCoord]:
    coords: dict[str, GeneCoord] = {}
    if not gff.exists():
        return coords
    for raw in gff.read_text().splitlines():
        if not raw or raw.startswith("#"):
            continue
        parts = raw.split("\t")
        if len(parts) < 9:
            continue
        seqid, _source, feature, start_raw, end_raw, _score, strand, _phase, attrs = parts
        if feature != "CDS":
            continue
        attr_map = _parse_gff_attrs(attrs)
        gene_id = attr_map.get("ID") or attr_map.get("Name")
        if not gene_id:
            continue
        try:
            start = int(start_raw)
            end = int(end_raw)
        except ValueError:
            continue
        coord = GeneCoord(contig=seqid, start=start, end=end, strand=strand)
        coords[gene_id] = coord
        prodigal_protein_id = _prodigal_protein_id(seqid=seqid, gene_id=gene_id)
        if prodigal_protein_id is not None:
            coords[prodigal_protein_id] = coord
    return coords


def _parse_gff_attrs(attrs: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in attrs.split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        out[key] = value
    return out


def _prodigal_protein_id(*, seqid: str, gene_id: str) -> Optional[str]:
    if "_" not in gene_id:
        return None
    _prefix, suffix = gene_id.split("_", 1)
    if not suffix:
        return None
    return f"{seqid}_{suffix}"


def _is_tigr_mtase_hit(hit: DomainHit) -> bool:
    text = _hit_text(hit)
    if _is_likely_rease_hit(hit):
        return False
    if _has_any(text, ("hsdm", "dam", "dcm", "type iii mod")):
        return True
    return _has_any(text, ("methyltransferase", "methylase", "mod subunit"))


def _is_methyltransferase_hit(hit: DomainHit) -> bool:
    if hit.source == "pfam" and hit.model_acc in MTASE_CATALYTIC_PFAMS:
        return True
    return _is_tigr_mtase_hit(hit)


def _is_likely_rease_hit(hit: DomainHit) -> bool:
    text = _hit_text(hit)
    if "methyltransferase" in text or "methylase" in text:
        return False
    if _has_any(text, ("hsdr", "type iii res", "res subunit")):
        return True
    if _has_any(text, ("restriction endonuclease", "restriction enzyme", "restriction subunit")):
        return True
    if "endonuclease" in text:
        return True
    return False


def _hit_text(hit: DomainHit) -> str:
    text = f"{hit.model_name} {hit.model_acc}".lower()
    text = text.replace("_", " ").replace("-", " ").replace("/", " ")
    return " ".join(text.split())


def _has_any(text: str, needles: tuple[str, ...]) -> bool:
    return any(needle in text for needle in needles)


def confidence_from_evalue(evalue: float) -> float:
    e = max(float(evalue), 1e-300)
    score = (-math.log10(e) - 5.0) / 50.0
    return max(0.0, min(1.0, score))


def _row_sort_key(row: dict[str, str]) -> tuple[str, int, str]:
    contig = row.get("contig") or "~"
    start_raw = row.get("start") or ""
    try:
        start = int(start_raw)
    except Exception:
        start = 1_000_000_000
    return (contig, start, row.get("candidate_id") or "")


def _safe_int(value: Optional[str]) -> int:
    try:
        return int(value or "0")
    except Exception:
        return 0
