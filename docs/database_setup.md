# Database Setup

`mtase-motif` ships the database management code, but it does not ship the
downloaded Pfam, TIGRFAMs, or REBASE payloads inside the package. Build the
runtime database after installation and keep it in a local cache directory.

## Quick Start

```bash
mtase-motif db init
mtase-motif db fetch pfam
mtase-motif db fetch rebase
mtase-motif db index
mtase-motif db status
```

Default DB directory: `~/.cache/mtase-motif/db`

## Offline Or Local Mirror Setup

Use `--source` when you already have local copies of the database payloads:

```bash
mtase-motif db fetch pfam --source /path/to/Pfam-A.hmm.gz
mtase-motif db fetch rebase --source /path/to/rebase_emboss_dir
```

Offline REBASE note: `db fetch rebase --source ...` only stages the normalized
REBASE tables. You still need `rebase_proteins.faa` or one or more REBASE
`*_Protein.txt` protein dumps in the source directory before `mtase-motif db
index` can build the REBASE protein search index.

Optional TIGRFAMs import is local-source only and always requires `--source`:

```bash
mtase-motif db fetch tigrfams --source /path/to/TIGRFAMs_15.0_HMM.LIB.gz
```

After any fetch or import, rebuild indexes:

```bash
mtase-motif db index
```

## What Gets Stored Locally

- `hmms/pfam/`: curated Pfam subset plus `hmmpress` outputs
- `hmms/tigrfams/`: optional TIGRFAM subset plus `hmmpress` outputs
- `rebase/<version>/`: normalized REBASE tables, protein FASTA, and local search indexes
- `manifests/`: provider metadata used by `db status`
