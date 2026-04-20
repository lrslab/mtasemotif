# `nanomotif` Comparison Notes

Date: 2026-03-29

Primary upstream source: [MicrobialDarkMatter/nanomotif](https://github.com/MicrobialDarkMatter/nanomotif)

## Executive Summary

`mtase-motif` and `nanomotif` are adjacent, but they are not solving the same problem in the same way.

- `mtase-motif` is sequence-first. The current pipeline discovers MTase candidates from protein/domain evidence, searches REBASE, optionally uses structure evidence, then scans the genome and writes QC outputs. The core path is in [`mtase_motif/runner.py`](/Users/li/Library/CloudStorage/Dropbox/rustproject/mtasemotif/mtase_motif/runner.py), [`mtase_motif/candidate_discovery.py`](/Users/li/Library/CloudStorage/Dropbox/rustproject/mtasemotif/mtase_motif/candidate_discovery.py), [`mtase_motif/motif_inference.py`](/Users/li/Library/CloudStorage/Dropbox/rustproject/mtasemotif/mtase_motif/motif_inference.py), and [`mtase_motif/scanning.py`](/Users/li/Library/CloudStorage/Dropbox/rustproject/mtasemotif/mtase_motif/scanning.py).
- `nanomotif` is methylation-evidence-first. Its main workflow starts from Nanopore/modkit pileups and discovers methylated motifs directly, then links those motifs back to MTases. The most relevant modules are `nanomotif/find_motifs_bin.py`, `nanomotif/motif.py`, `nanomotif/postprocess.py`, and `nanomotif/mtase_linker/`.

The closest overlap is not the whole `nanomotif` package. It is specifically their `MTase-linker` subworkflow plus a few reusable ideas from their motif discovery engine.

## Where `nanomotif` Looks Strong

### 1. Methylation-backed motif support

`nanomotif` carries actual methylated vs non-methylated counts through its outputs (`n_mod`, `n_nomod`, `mod_type`, `mod_position`, motif complement). That is stronger evidence than our current sequence transfer plus FIMO/QC-only approach.

Implication for us:

- add an optional pileup-backed validation stage after motif inference
- use observed methylation evidence to confirm or downgrade transferred motifs
- expose motif support counts in final outputs

This is the highest-value feature to borrow.

### 2. Probabilistic scoring instead of threshold-only confidence

`nanomotif/model.py` uses a Beta-Bernoulli model, and `nanomotif/find_motifs_bin.py` uses posterior-predictive style scoring when refining motif candidates and deciding whether submotifs should survive cleanup.

Implication for us:

- replace some of the current hard-threshold confidence logic with an evidence score
- combine REBASE homology, related-candidate transfer, structure hints, and methylation support into a more defensible confidence model

We do not need their exact de novo search to benefit from this.

### 3. Motif post-processing and normalization

`nanomotif/postprocess.py` and `nanomotif/motif.py` do useful cleanup:

- remove noisy motifs
- collapse reverse complements
- reason about parent/child submotif relationships
- keep motif sequence and modified-base position together

Implication for us:

- add a post-inference deduplication layer in [`mtase_motif/motif_inference.py`](/Users/li/Library/CloudStorage/Dropbox/rustproject/mtasemotif/mtase_motif/motif_inference.py)
- reduce duplicate or near-duplicate motif calls from REBASE transfer, related-candidate rescue, and structure-panel rescue
- track modified-base position explicitly in final motif outputs

### 4. Better ambiguity reporting

`nanomotif`’s `MTase-linker` reports both confident links and unresolved candidate genes. That is better than forcing everything into a single best assignment.

Implication for us:

- add first-class ambiguous-link reporting
- add a dedicated assignment-style output, for example `motif_assignment.tsv`, instead of burying ambiguity in secondary fields
- preserve multiple plausible motif or gene assignments when evidence is insufficient to collapse to one answer
- make the summary output clearer about “best call” vs “candidate set”

This matches the existing shape of our inference code, which already carries related candidates and hint fields, but does not surface ambiguity as cleanly as it could.

### 5. A stronger motif result abstraction

`nanomotif/motif.py` uses `Motif`, `MotifTree`, and `MotifSearchResult` to centralize motif representation, derived fields, and output formatting.

Implication for us:

- consider a lightweight local result model for motif calls and summaries
- avoid scattering TSV schema assumptions across inference, scanning, and QC code
- make downstream additions like complement fields or support counts less fragile

### 6. Locus-aware assignment heuristics

One useful idea from `MTase-linker` is that linking is not purely motif-to-gene similarity. It also uses RM-system context to reduce ambiguity.

Implication for us:

- promote local neighborhood context from [`mtase_motif/candidate_discovery.py`](/Users/li/Library/CloudStorage/Dropbox/rustproject/mtasemotif/mtase_motif/candidate_discovery.py) into an explicit reusable field such as `locus_group_id` or `system_id`
- let downstream inference prefer assignments that are consistent within a local RM locus
- separate “good motif hit” from “good motif hit in the right local system context”

### 7. Motif-architecture agreement as a confidence feature

`nanomotif` explicitly reasons about motif architecture such as palindromic, non-palindromic, and bipartite calls.

Implication for us:

- score agreement between predicted MTase type/subtype and motif architecture
- use architecture agreement as an explicit confidence feature, not just a warning in the final summary
- downgrade assignments that are sequence-similar but architecture-inconsistent

## What Fits Poorly

### 1. Full de novo motif discovery as the default workflow

`nanomotif`’s core search depends on:

- Nanopore/modkit pileups
- methylation frequency thresholds
- bin-aware processing
- metagenomic context

That is not the current design center of `mtase-motif`, which is a single-genome, local-first, sequence-first MTase workflow.

Conclusion:

- do not replace the current core pipeline with `nanomotif`-style de novo search
- if implemented, make de novo methylation support an optional evidence module

### 2. Bin contamination and unbinned-contig inclusion

The `binnary/` subpackage is useful for methylation-guided metagenomic cleanup, but it is outside our current scope.

Conclusion:

- do not import `detect_contamination` or `include_contigs` into `mtase-motif`

### 3. The `MTase-linker` implementation style

The `MTase-linker` ideas are relevant, but the implementation relies on:

- extra Conda environments
- Snakemake orchestration
- DefenseFinder
- pandas scripts with fairly ad hoc rule logic

Conclusion:

- borrow heuristics, not the implementation
- keep our current Python-first package structure

## What `mtase-motif` Already Does Well

There are areas where `mtase-motif` already has a cleaner or more directly usable design:

- local-first DB management and provider indexing in [`mtase_motif/db/`](/Users/li/Library/CloudStorage/Dropbox/rustproject/mtasemotif/mtase_motif/db)
- a direct single-genome run path in [`mtase_motif/runner.py`](/Users/li/Library/CloudStorage/Dropbox/rustproject/mtasemotif/mtase_motif/runner.py)
- optional structure-assisted rescue and template-panel inference in [`mtase_motif/structure_motif/template_panel.py`](/Users/li/Library/CloudStorage/Dropbox/rustproject/mtasemotif/mtase_motif/structure_motif/template_panel.py)
- exact-IUPAC fallback and explicit QC output in [`mtase_motif/qc.py`](/Users/li/Library/CloudStorage/Dropbox/rustproject/mtasemotif/mtase_motif/qc.py)

This suggests we should treat `nanomotif` as a source of complementary evidence ideas, not as a package to mirror.

## Recommended Incorporation Plan

### Priority 1

Add optional methylation evidence input and validation.

Concrete direction:

- accept modkit pileup or a normalized methylation-support table as optional input
- validate inferred motifs against observed methylated and non-methylated motif instances
- add support counts and evidence-based confidence updates to motif outputs

### Priority 2

Improve motif output semantics.

Concrete direction:

- add modified-base position to final motif outputs
- add motif class fields such as palindromic / non-palindromic / bipartite
- expose complement or reverse-complement normalization where relevant
- report ambiguous candidate sets separately from the primary call
- add a dedicated assignment table with explicit `linked`, `ambiguous`, and `unresolved` states

### Priority 3

Add post-processing for related motif calls.

Concrete direction:

- collapse complements
- remove redundant submotifs where evidence does not support keeping both
- merge duplicate calls arising from multiple inference routes
- resolve competing motif-to-gene assignments after the initial transfer step instead of only during candidate scoring

### Priority 4

Introduce locus-aware assignment and a lightweight motif result abstraction.

Concrete direction:

- promote neighborhood-derived grouping into a stable field such as `locus_group_id`
- use that grouping during ambiguity resolution and summary reporting
- define a shared result type or schema object for motif calls
- centralize derived-field generation and writing logic
- reduce duplication between inference, scanning, and summary generation

### Priority 5

Only after the above, evaluate a limited de novo discovery mode.

Concrete direction:

- keep it optional
- scope it to validation or rescue for unresolved candidates
- avoid dragging metagenomic bin logic into the default workflow

### Priority 6

Add workflow diagnostics for environment and evidence quality.

Concrete direction:

- add a lightweight `mtase-motif doctor` or equivalent diagnostic command
- report missing databases, missing executables, weak evidence layers, and contradictory motif-assignment signals in one place
- make it easier to tell whether a weak result is caused by biology, missing inputs, or runtime setup

## Suggested Near-Term Tasks

1. Extend motif outputs with `mod_position`, motif class, and explicit `linked` vs `ambiguous` reporting.
2. Add a dedicated assignment table that separates primary calls from alternate candidates.
3. Promote neighborhood context into a reusable `locus_group_id` or equivalent field.
4. Add a post-processing layer to collapse complementary and redundant motif calls.
5. Design a minimal methylation-support input format that can be consumed without importing `nanomotif`.
6. Prototype a simple evidence score that combines homology, structure, architecture agreement, and methylation support.

## Bottom Line

The main thing worth borrowing from `nanomotif` is not its full pipeline. It is the idea that motif calls should be backed by observed methylation evidence and normalized with stronger motif-aware post-processing.

For `mtase-motif`, the best path is:

- keep the current sequence-first architecture
- add optional methylation-backed validation
- improve result representation and ambiguity reporting
- avoid taking on metagenomic bin workflows unless the package scope changes
