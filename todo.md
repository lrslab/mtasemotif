# TODO

## Review Findings (2026-03-29)

Reference: correctness/docs/usage review across the current package.
Validation note: `pytest -q` passed locally (`79 passed`), so the items below are mostly edge cases, optional paths, and cleanup gaps that the current tests do not cover.

### High priority

- [ ] Fix structure-file matching so candidate IDs cannot collide by prefix during fallback lookup (for example `gene1` resolving to `gene10_model.pdb`).
- [ ] Fix template-panel provenance and confidence reporting so the chosen `best_target` reflects the hit(s) that actually contributed DNA to the panel.
- [ ] Reject `--foldseek-labels` when no Foldseek DB is available, or add a real execution path that can consume labels without a DB-backed Foldseek search.
- [ ] Fix CLI top-level error handling so user-facing failures exit cleanly without a Python traceback and with the intended non-zero status code.

### Medium priority

- [ ] Fix structure-training readiness accounting so `usable_for_template_panel` only becomes true when panel construction actually succeeds.
- [ ] Make the offline REBASE docs accurate: `db fetch rebase --source ...` is not enough on its own unless `rebase_proteins.faa` is also provided for indexing.
- [ ] Make `db fetch` help text explicit that `tigrfams` is local-source-only and that `--source` is required for that target.
- [ ] Document the Foldseek/structure-mode prerequisites more clearly in the CLI help and README, including the `--structures-dir` precondition for Foldseek-related options.
- [ ] Make the end-to-end README walkthrough self-contained or clearly mark the example genome path as machine-local rather than checkout-local.

### Cleanup / Dead Code

- [ ] Remove or wire up `merge_rows()` in `mtase_motif/pipeline_support.py`.
- [ ] Remove or wire up `layout()` in `mtase_motif/db/layout.py`.
- [ ] Remove or wire up `PanelResult` and `build_weighted_kmer_panel()` in `mtase_motif/panel.py`.
- [ ] Either use or drop the `Provider.status()` / `Provider.validate()` interface in the current CLI flow.
- [ ] Remove the unused `platformdirs` dependency from packaging metadata if it is no longer needed.

Possible improvements identified from the `nanomotif` comparison.
Reference: [docs/nanomotif_comparison.md](/Users/li/Library/CloudStorage/Dropbox/rustproject/mtasemotif/docs/nanomotif_comparison.md)

## High priority

- [ ] Add optional methylation-evidence input and validation.
- [ ] Accept `modkit` pileup or a normalized methylation-support table as optional input.
- [ ] Validate inferred motifs against observed methylated and non-methylated motif instances.
- [ ] Add support counts and evidence-based confidence updates to motif outputs.

- [ ] Improve motif output semantics.
- [ ] Add modified-base position (`mod_position`) to final motif outputs.
- [ ] Add motif class fields such as palindromic, non-palindromic, and bipartite.
- [ ] Normalize complements or reverse complements where relevant.
- [ ] Add explicit `linked`, `ambiguous`, and `unresolved` assignment states.
- [ ] Add a dedicated `motif_assignment.tsv` that separates primary calls from alternate candidates.

## Medium priority

- [ ] Add post-processing for related motif calls.
- [ ] Collapse complementary motif calls.
- [ ] Remove redundant submotifs where evidence does not support keeping both.
- [ ] Merge duplicate calls arising from REBASE transfer, related-candidate rescue, and structure-panel rescue.
- [ ] Add a post-inference ambiguity-resolution stage for competing motif-to-gene assignments.

- [ ] Make assignment more locus-aware.
- [ ] Promote neighborhood-derived grouping into a stable `locus_group_id` or `system_id`.
- [ ] Use local RM-locus context during ambiguity resolution and summary reporting.
- [ ] Separate “good motif hit” from “good motif hit in the right local system context”.

- [ ] Make motif-architecture agreement an explicit confidence feature.
- [ ] Score agreement between predicted MTase type/subtype and motif architecture.
- [ ] Downgrade assignments that are sequence-similar but architecture-inconsistent.

## Lower priority

- [ ] Introduce a lightweight shared result abstraction for motif calls and summaries.
- [ ] Centralize derived-field generation and output formatting instead of spreading TSV schema logic across modules.
- [ ] Add complement fields and support-count fields through a single result schema.

- [ ] Evaluate a limited de novo discovery mode as an optional rescue path.
- [ ] Keep any de novo methylation workflow optional, not part of the default sequence-first pipeline.
- [ ] Scope it to validation or rescue for unresolved candidates.

- [ ] Add a lightweight `mtase-motif doctor` command.
- [ ] Report missing executables, missing databases, weak evidence layers, and contradictory assignment signals in one place.

## Non-goals for now

- [ ] Do not import `nanomotif` metagenomic bin contamination workflows.
- [ ] Do not import `nanomotif` unbinned-contig inclusion workflows.
- [ ] Do not copy the separate Snakemake/Conda dependency island from `MTase-linker`.
