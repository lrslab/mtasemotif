from __future__ import annotations

import logging
from dataclasses import asdict
from enum import Enum
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from mtase_motif import __version__
from mtase_motif.config import DEFAULT_DB_DIR, DbLayout, RunConfig, expand_path
from mtase_motif.db.manifest import (
    DbManifest,
    load_manifest,
    load_provider_manifest,
    write_provider_manifest,
    write_manifest,
)
from mtase_motif.db.providers import get_provider_class, list_providers
from mtase_motif.logging_utils import configure_logging
from mtase_motif.runner import build_run_plan, resolve_run_context, run_pipeline
from mtase_motif.util import MtaseMotifError, which


app = typer.Typer(add_completion=False)
db_app = typer.Typer(add_completion=False, no_args_is_help=True)
app.add_typer(db_app, name="db")

console = Console()
log = logging.getLogger("mtase_motif")


class FetchTarget(str, Enum):
    pfam = "pfam"
    tigrfams = "tigrfams"
    rebase = "rebase"


def validate_runtime_dependencies(cfg: RunConfig, *, dry_run: bool) -> None:
    if cfg.foldseek_db is not None and cfg.structures_dir is None:
        raise MtaseMotifError("--foldseek-db requires --structures-dir")
    default_foldseek_db = cfg.db_dir / "structures" / "pdb" / "foldseek_db"
    available_foldseek_db = cfg.foldseek_db
    if available_foldseek_db is None and default_foldseek_db.exists():
        available_foldseek_db = default_foldseek_db
    if cfg.foldseek_labels is not None:
        if cfg.structures_dir is None:
            raise MtaseMotifError("--foldseek-labels requires --structures-dir")
        if available_foldseek_db is None:
            raise MtaseMotifError(
                "--foldseek-labels requires a Foldseek DB from --foldseek-db or the "
                f"default path {default_foldseek_db}; label transfer only runs when "
                "Foldseek search is enabled."
            )

    if dry_run:
        return

    if cfg.proteins_faa is None and which("prodigal") is None:
        raise MtaseMotifError("prodigal not found in PATH (required unless --proteins is provided)")
    if which("hmmscan") is None:
        raise MtaseMotifError("hmmscan not found in PATH (HMMER required for candidate detection)")
    tigr_subset = cfg.db_dir / "hmms" / "tigrfams" / "subset.hmm"
    if tigr_subset.exists() and which("hmmsearch") is None:
        raise MtaseMotifError(
            f"Found TIGRFAMs subset but hmmsearch is not in PATH: {tigr_subset}\n"
            "Install HMMER tools including hmmsearch, or remove the TIGRFAMs subset."
        )
    foldseek_required = cfg.structures_dir is not None and available_foldseek_db is not None
    if foldseek_required and which("foldseek") is None:
        raise MtaseMotifError(
            "foldseek not found in PATH (required when using --foldseek-db or the "
            f"auto-discovered default DB at {default_foldseek_db})"
        )
    if which("fimo") is None:
        raise MtaseMotifError("fimo not found in PATH (MEME suite required for motif scanning)")


def validate_pressed_hmm_db(hmm_path: Path, *, label: str) -> None:
    missing = [Path(str(hmm_path) + ext) for ext in (".h3f", ".h3i", ".h3m", ".h3p")]
    missing = [p for p in missing if not p.exists()]
    if missing:
        raise MtaseMotifError(
            f"{label} is missing hmmpress index files for {hmm_path}.\n"
            "Run: mtase-motif db index"
        )


@app.callback(invoke_without_command=True)
def main_callback(
    ctx: typer.Context,
    version: bool = typer.Option(False, "--version", help="Show version and exit."),
    log_level: str = typer.Option("INFO", "--log-level", help="Logging level."),
) -> None:
    if version:
        console.print(__version__)
        raise typer.Exit(code=0)
    ctx.ensure_object(dict)
    ctx.obj["log_level"] = log_level
    configure_logging(level=log_level)
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise typer.Exit(code=0)


@db_app.command("init")
def db_init(
    db_dir: Path = typer.Option(DEFAULT_DB_DIR, "--db-dir", help="Database directory."),
    force: bool = typer.Option(False, "--force", help="Overwrite existing manifest."),
) -> None:
    db_dir = expand_path(db_dir)
    layout = DbLayout(db_dir=db_dir)
    layout.ensure_dirs()

    if layout.manifest_path.exists() and not force:
        raise MtaseMotifError(
            f"DB already initialized (manifest exists): {layout.manifest_path}\n"
            "Re-run with --force to overwrite."
        )

    write_manifest(layout.manifest_path, DbManifest())
    console.print(f"Initialized db at {db_dir}")


@db_app.command("status")
def db_status(
    db_dir: Path = typer.Option(DEFAULT_DB_DIR, "--db-dir", help="Database directory."),
) -> None:
    db_dir = expand_path(db_dir)
    layout = DbLayout(db_dir=db_dir)
    if not layout.manifest_path.exists():
        raise MtaseMotifError(
            f"DB not initialized: {layout.manifest_path}\nRun: mtase-motif db init"
        )
    manifest = load_manifest(layout.manifest_path)
    providers: dict[str, object] = {}
    for provider_name in list_providers():
        provider = get_provider_class(provider_name)(layout)
        info = dict(provider.status())
        try:
            provider.validate()
            info["valid"] = True
        except MtaseMotifError as e:
            info["valid"] = False
            info["validation_error"] = str(e)

        manifest_path = layout.manifests_dir / f"{provider_name}.json"
        if manifest_path.exists():
            try:
                info["manifest"] = asdict(load_provider_manifest(manifest_path))
            except Exception as e:
                info["manifest_error"] = str(e)
                info["manifest_path"] = str(manifest_path)
        providers[provider_name] = info
    console.print_json(data={"manifest": asdict(manifest), "providers": providers})


@db_app.command("fetch")
def db_fetch(
    target: FetchTarget = typer.Argument(
        ...,
        help="Which database to fetch. tigrfams is local-source-only and requires --source.",
    ),
    db_dir: Path = typer.Option(DEFAULT_DB_DIR, "--db-dir", help="Database directory."),
    source: Optional[Path] = typer.Option(
        None, "--source", help="Optional local source path (required for tigrfams)."
    ),
) -> None:
    db_dir = expand_path(db_dir)
    layout = DbLayout(db_dir=db_dir)
    if not layout.manifest_path.exists():
        raise MtaseMotifError("DB not initialized. Run: mtase-motif db init")

    if target is FetchTarget.tigrfams and source is None:
        raise MtaseMotifError(
            "tigrfams fetch is local-source-only and requires --source.\n"
            "Example: mtase-motif db fetch tigrfams --source /path/to/TIGRFAMs_15.0_HMM.LIB.gz"
        )

    if source is not None:
        source = expand_path(source)

    layout.manifests_dir.mkdir(parents=True, exist_ok=True)
    provider_cls = get_provider_class(target.value)
    provider = provider_cls(layout)
    pm = provider.fetch(source=source)

    provider_manifest_path = layout.manifests_dir / f"{target.value}.json"
    write_provider_manifest(provider_manifest_path, pm)

    manifest = load_manifest(layout.manifest_path)
    rel = str(provider_manifest_path.relative_to(layout.db_dir))
    manifest.providers[target.value] = {"manifest": rel}
    write_manifest(layout.manifest_path, manifest)

    console.print(f"Fetched {target.value} into {db_dir}")


@db_app.command("index")
def db_index(
    db_dir: Path = typer.Option(DEFAULT_DB_DIR, "--db-dir", help="Database directory."),
) -> None:
    db_dir = expand_path(db_dir)
    layout = DbLayout(db_dir=db_dir)
    if not layout.manifest_path.exists():
        raise MtaseMotifError("DB not initialized. Run: mtase-motif db init")

    manifest = load_manifest(layout.manifest_path)
    if not layout.manifests_dir.exists():
        console.print("No provider manifests found to index.")
        raise typer.Exit(code=0)

    errors: list[tuple[str, str]] = []
    for path in sorted(layout.manifests_dir.glob("*.json")):
        provider_name = path.stem
        provider_cls = get_provider_class(provider_name)
        provider = provider_cls(layout)
        try:
            pm = provider.index()
            write_provider_manifest(path, pm)
            manifest.providers.setdefault(provider_name, {"manifest": str(path.relative_to(layout.db_dir))})
            console.print(f"Indexed {provider_name}")
        except MtaseMotifError as e:
            errors.append((provider_name, str(e)))
            console.print(f"[red]Index failed[/red] {provider_name}: {e}")

    write_manifest(layout.manifest_path, manifest)
    if errors:
        raise MtaseMotifError(
            "Some providers failed to index:\n"
            + "\n".join([f"- {name}: {msg}" for name, msg in errors])
        )
    console.print(f"Indexed providers in {db_dir}")


@app.command("run")
def run(
    ctx: typer.Context,
    genome: Path = typer.Option(..., "--genome", exists=True, help="Genome FASTA (.fna/.fa)."),
    out: Path = typer.Option(Path("results"), "--out", help="Output directory."),
    db_dir: Path = typer.Option(DEFAULT_DB_DIR, "--db-dir", help="Database directory."),
    proteins: Optional[Path] = typer.Option(
        None, "--proteins", exists=True, help="Optional proteins FASTA; skips gene calling."
    ),
    structures_dir: Optional[Path] = typer.Option(
        None,
        "--structures-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        help="Optional directory of candidate structures (PDB/mmCIF). Required for Foldseek and structure-template mode.",
    ),
    foldseek_db: Optional[Path] = typer.Option(
        None,
        "--foldseek-db",
        exists=True,
        help=(
            "Optional Foldseek target DB prefix. Requires --structures-dir when set; "
            "otherwise the CLI auto-discovers <db-dir>/structures/pdb/foldseek_db when present."
        ),
    ),
    foldseek_labels: Optional[Path] = typer.Option(
        None,
        "--foldseek-labels",
        exists=True,
        help=(
            "Optional TSV mapping Foldseek target_id -> motif_iupac. Requires --structures-dir "
            "and an explicit or auto-discovered Foldseek DB."
        ),
    ),
    mmcif_mirror_dir: Optional[Path] = typer.Option(
        None,
        "--mmcif-mirror-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        help="Optional directory containing PDB mmCIF files for template extraction.",
    ),
    mmcif_cache_dir: Optional[Path] = typer.Option(
        None,
        "--mmcif-cache-dir",
        file_okay=False,
        dir_okay=True,
        help="Directory to cache downloaded mmCIF templates (default: <out>/work/structures/mmcif).",
    ),
    mmcif_download: bool = typer.Option(
        False,
        "--download-mmcif/--no-download-mmcif",
        help="Allow downloading missing PDB mmCIF templates during structure-panel inference.",
    ),
    template_top_hits: int = typer.Option(
        3,
        "--template-top-hits",
        min=1,
        help="How many top Foldseek hits to consider for template-panel inference.",
    ),
    template_min_alntmscore: float = typer.Option(
        0.5,
        "--template-min-alntmscore",
        min=0.0,
        max=1.0,
        help="Minimum Foldseek alignment TM-score for template-panel inference.",
    ),
    panel_kmin: int = typer.Option(
        6, "--panel-kmin", min=4, help="Minimum k-mer length for template-panel inference."
    ),
    panel_kmax: int = typer.Option(
        8, "--panel-kmax", min=4, help="Maximum k-mer length for template-panel inference."
    ),
    panel_size: int = typer.Option(
        25, "--panel-size", min=1, help="Number of k-mers to use when building panel PWM."
    ),
    jobs: int = typer.Option(1, "-j", "--jobs", min=1, help="Parallel jobs for native tools."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print the native pipeline plan and exit."),
) -> None:
    cfg = RunConfig(
        genome_fasta=genome,
        db_dir=db_dir,
        out_dir=out,
        proteins_faa=proteins,
        structures_dir=structures_dir,
        foldseek_db=foldseek_db,
        foldseek_labels=foldseek_labels,
        mmcif_mirror_dir=mmcif_mirror_dir,
        mmcif_cache_dir=mmcif_cache_dir,
        mmcif_download=mmcif_download,
        template_top_hits=template_top_hits,
        template_min_alntmscore=template_min_alntmscore,
        panel_kmin=panel_kmin,
        panel_kmax=panel_kmax,
        panel_size=panel_size,
    )
    cfg = cfg.normalized()
    log_level = "INFO"
    if ctx.obj and isinstance(ctx.obj, dict):
        log_level = str(ctx.obj.get("log_level", log_level))

    layout = DbLayout(db_dir=cfg.db_dir)
    if not layout.manifest_path.exists():
        raise MtaseMotifError(f"DB not initialized: {layout.manifest_path}\nRun: mtase-motif db init")

    validate_runtime_dependencies(cfg, dry_run=dry_run)
    run_ctx = resolve_run_context(cfg)

    if not run_ctx.pfam_subset.exists():
        raise MtaseMotifError(
            f"Missing Pfam subset HMM: {run_ctx.pfam_subset}\n"
            "Run: mtase-motif db fetch pfam && mtase-motif db index\n"
            "(or use --source <Pfam-A.hmm(.gz)> for offline mode)"
        )
    if not dry_run:
        validate_pressed_hmm_db(run_ctx.pfam_subset, label="Pfam subset HMM")

    if not run_ctx.rebase_enzymes.exists():
        raise MtaseMotifError(
            f"Missing {run_ctx.rebase_enzymes}. Re-run: mtase-motif db fetch rebase ..."
        )
    if not run_ctx.rebase_proteins.exists():
        raise MtaseMotifError(
            f"Missing {run_ctx.rebase_proteins}.\n"
            "Re-run: mtase-motif db fetch rebase && mtase-motif db index\n"
            "(or place rebase_proteins.faa or one or more REBASE *_Protein.txt protein "
            "dumps in your rebase --source dir and re-run fetch/index)."
        )

    if not dry_run:
        if run_ctx.tigr_subset is not None:
            validate_pressed_hmm_db(run_ctx.tigr_subset, label="TIGRFAMs subset HMM")
        mmseqs_ok = which("mmseqs") is not None and (
            run_ctx.mmseqs_db.exists() or Path(str(run_ctx.mmseqs_db) + ".dbtype").exists()
        )
        blast_ok = which("blastp") is not None and any(
            Path(str(run_ctx.blast_db) + ext).exists() for ext in (".pin", ".psq", ".phr")
        )
        if not (mmseqs_ok or blast_ok):
            raise MtaseMotifError(
                f"No usable REBASE protein index found under {run_ctx.rebase_version_dir}.\n"
                "Run: mtase-motif db index (requires mmseqs or BLAST+ tools)."
            )

    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    configure_logging(level=log_level, log_file=cfg.out_dir / "logs" / "run.jsonl")
    log.info(
        "run_config genome=%s db_dir=%s out_dir=%s proteins=%s",
        cfg.genome_fasta,
        cfg.db_dir,
        cfg.out_dir,
        cfg.proteins_faa,
    )

    if dry_run:
        console.print("Planned native pipeline:")
        for idx, step in enumerate(build_run_plan(run_ctx), start=1):
            console.print(f"{idx}. {step}")
        raise typer.Exit(code=0)

    run_pipeline(run_ctx, jobs=jobs, dry_run=False, logger=log)


def main() -> None:
    try:
        app()
    except MtaseMotifError as e:
        console.print(f"[red]Error:[/red] {e}")
        raise SystemExit(2)
    except typer.Exit:
        raise
    except Exception as e:
        console.print(f"[red]Unexpected error:[/red] {e}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
