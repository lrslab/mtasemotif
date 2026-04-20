from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from mtase_motif.util import MtaseMotifError


DEFAULT_DB_DIR = Path.home() / ".cache" / "mtase-motif" / "db"


def expand_path(path: Path) -> Path:
    return Path(path).expanduser().resolve()


@dataclass(frozen=True)
class DbLayout:
    db_dir: Path

    @property
    def manifest_path(self) -> Path:
        return self.db_dir / "manifest.json"

    @property
    def manifests_dir(self) -> Path:
        return self.db_dir / "manifests"

    @property
    def hmms_dir(self) -> Path:
        return self.db_dir / "hmms"

    @property
    def pfam_dir(self) -> Path:
        return self.hmms_dir / "pfam"

    @property
    def tigrfams_dir(self) -> Path:
        return self.hmms_dir / "tigrfams"

    @property
    def rebase_dir(self) -> Path:
        return self.db_dir / "rebase"

    @property
    def structures_dir(self) -> Path:
        return self.db_dir / "structures"

    @property
    def pdb_structures_dir(self) -> Path:
        return self.structures_dir / "pdb"

    def ensure_dirs(self) -> None:
        self.db_dir.mkdir(parents=True, exist_ok=True)
        self.manifests_dir.mkdir(parents=True, exist_ok=True)
        self.pfam_dir.mkdir(parents=True, exist_ok=True)
        self.tigrfams_dir.mkdir(parents=True, exist_ok=True)
        self.rebase_dir.mkdir(parents=True, exist_ok=True)
        self.pdb_structures_dir.mkdir(parents=True, exist_ok=True)


@dataclass(frozen=True)
class RunConfig:
    genome_fasta: Path
    db_dir: Path
    out_dir: Path
    proteins_faa: Optional[Path] = None
    structures_dir: Optional[Path] = None
    foldseek_db: Optional[Path] = None
    foldseek_labels: Optional[Path] = None
    mmcif_mirror_dir: Optional[Path] = None
    mmcif_cache_dir: Optional[Path] = None
    mmcif_download: bool = False
    template_top_hits: int = 3
    template_min_alntmscore: float = 0.5
    panel_kmin: int = 6
    panel_kmax: int = 8
    panel_size: int = 25

    def normalized(self) -> "RunConfig":
        panel_kmin = int(self.panel_kmin)
        panel_kmax = int(self.panel_kmax)
        if panel_kmin > panel_kmax:
            raise MtaseMotifError("--panel-kmin must be less than or equal to --panel-kmax")

        return RunConfig(
            genome_fasta=expand_path(self.genome_fasta),
            db_dir=expand_path(self.db_dir),
            out_dir=expand_path(self.out_dir),
            proteins_faa=expand_path(self.proteins_faa) if self.proteins_faa else None,
            structures_dir=expand_path(self.structures_dir) if self.structures_dir else None,
            foldseek_db=expand_path(self.foldseek_db) if self.foldseek_db else None,
            foldseek_labels=expand_path(self.foldseek_labels) if self.foldseek_labels else None,
            mmcif_mirror_dir=expand_path(self.mmcif_mirror_dir) if self.mmcif_mirror_dir else None,
            mmcif_cache_dir=expand_path(self.mmcif_cache_dir) if self.mmcif_cache_dir else None,
            mmcif_download=bool(self.mmcif_download),
            template_top_hits=int(self.template_top_hits),
            template_min_alntmscore=float(self.template_min_alntmscore),
            panel_kmin=panel_kmin,
            panel_kmax=panel_kmax,
            panel_size=int(self.panel_size),
        )
