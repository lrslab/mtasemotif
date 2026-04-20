# mtase-motif

`mtase-motif` is a local-first Python package for finding bacterial DNA
methyltransferase candidates and transferring or inferring methylation motifs
from a single genome.

The package keeps the database management code in `mtase_motif/`, but it does
not bundle downloaded Pfam, TIGRFAMs, or REBASE payloads in the published
artifact.

## Install

Install the Python package in a local virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev]'
```

For a PyPI release, install the distribution as:

```bash
python -m pip install mtasemotif
```

The installed CLI command remains:

```bash
mtase-motif --help
```

Install native executables separately. The repo includes an optional
`environment.yml` for a conda env that provides the non-Python tools:

```bash
conda env create -f environment.yml
conda activate mtase
```

## End-To-End Walkthrough

The commands below cover the full sequence-first workflow from a clean checkout
to a completed run on the example E. coli genome stored on this machine at
`/Users/li/data/hammerhead_test/ecoli.fa`. That path is machine-local, not part
of the checkout.

### 1. Prepare the environments

Activate the native-tool conda environment first, then the local Python venv:

```bash
conda activate mtase
source .venv/bin/activate
```

This keeps `prodigal`, `hmmscan`, `hmmsearch`, `hmmpress`, `mmseqs` or
`blastp`, and `fimo` on `PATH` while the Python package stays in `.venv`.
`tigrfams` remains local-source-only, so you only use `--source` for that
target when you already have a local TIGRFAMs HMM install.

### 2. Create and populate the local database

The walkthrough below uses a project-local DB directory so the full run is
self-contained:

```bash
DB_DIR=$PWD/db

mtase-motif db init --db-dir "$DB_DIR"
mtase-motif db fetch pfam --db-dir "$DB_DIR"
mtase-motif db fetch rebase --db-dir "$DB_DIR"
mtase-motif db index --db-dir "$DB_DIR"
mtase-motif db status --db-dir "$DB_DIR"
```

Downloaded databases live outside the package by default under
`~/.cache/mtase-motif/db`, but using `--db-dir "$DB_DIR"` keeps this walkthrough
local to the repo.

Offline and local-mirror examples:

```bash
mtase-motif db fetch pfam --db-dir "$DB_DIR" --source /path/to/Pfam-A.hmm.gz
mtase-motif db fetch rebase --db-dir "$DB_DIR" --source /path/to/rebase_emboss_dir
mtase-motif db fetch tigrfams --db-dir "$DB_DIR" --source /path/to/TIGRFAMs_15.0_HMM.LIB.gz
mtase-motif db index --db-dir "$DB_DIR"
```

For offline REBASE, `--source` alone is not enough. The source directory must
also contain `rebase_proteins.faa` or one or more REBASE `*_Protein.txt`
protein dumps before `mtase-motif db index` can build the REBASE protein
search index.

More database download and import notes are in `docs/database_setup.md`.

### 3. Run motif prediction on a genome FASTA

Use the example genome from `/Users/li/data/hammerhead_test`:

```bash
GENOME=/Users/li/data/hammerhead_test/ecoli.fa
OUT_DIR=$PWD/results/ecoli

mtase-motif run --genome "$GENOME" --db-dir "$DB_DIR" --out "$OUT_DIR" -j 4
```

Replace `"$GENOME"` with your own `.fa` or `.fna` file when running on a new
genome.

Optional structure-assisted runs require local candidate structures via
`--structures-dir`. Foldseek-backed steps also require `foldseek` on `PATH` and
either `--foldseek-db` or the default DB at
`<db-dir>/structures/pdb/foldseek_db`; `--foldseek-labels` only applies when
that Foldseek DB-backed search path is available.

If you already have predicted proteins, you can skip gene calling:

```bash
mtase-motif run --genome "$GENOME" --proteins /path/to/proteins.faa --db-dir "$DB_DIR" --out "$OUT_DIR" -j 4
```

### 4. Inspect the outputs

Core outputs in `"$OUT_DIR"`:

- `"$OUT_DIR"/mtase_candidates.tsv`
- `"$OUT_DIR"/motif_calls.tsv`
- `"$OUT_DIR"/motif_assignment.tsv`
- `"$OUT_DIR"/summary.tsv`
- `"$OUT_DIR"/<candidate_id>/motif/pwm.meme`
- `"$OUT_DIR"/<candidate_id>/fimo/fimo.tsv`
- `"$OUT_DIR"/<candidate_id>/qc/qc.json`

`motif_calls.tsv` now carries derived motif semantics such as `mod_position`,
`motif_class`, canonical/reverse-complement forms, and an overall
`assignment_state`. `motif_assignment.tsv` separates the primary call from
alternate related-candidate or hint routes when a candidate is unresolved or
ambiguous.

Quick checks:

```bash
ls "$OUT_DIR"
head "$OUT_DIR/mtase_candidates.tsv"
head "$OUT_DIR/summary.tsv"
```

## Runtime Tools

These executables must be on `PATH` for the sequence-first workflow:

- `prodigal`
- `hmmscan`, `hmmsearch`, `hmmpress`
- `mmseqs` or `blastp` plus `makeblastdb`
- `fimo`

## Development

```bash
make lint
make test
make sdist-check
make package-check
```

Tagged releases are set up for GitHub Releases and PyPI publishing through the
GitHub Actions workflows in `.github/workflows/`.
