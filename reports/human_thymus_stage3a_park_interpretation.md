# Park Human TEC Stage 3A Interpretation

Search date: 2026-07-06

## Local Parsing Status

- Park TEC H5AD parsed: yes.
- Dataset shape: 18,524 cells x 32,088 features.
- Raw object exists: yes.
- Layers found: none.
- Matrix/layer semantics: partial. `X` and `raw.X` are available as sparse float32 matrices, but this Stage 3A audit does not assert count or normalization semantics.

## Usable Metadata Fields

- Usable donor/sample field: `donor_id` for donor grouping, with `sample_id` retained as a sample/library-level field.
- Usable age/development-stage field: `development_stage`; `development_stage_ontology_term_id` is also available.
- Usable sex field: `sex`.
- Usable cell-type/subtype field: `celltypes`; `cell_type` is also available as a broader ontology-aligned label.

## LOX-Family Gene Presence

All target genes were detected by exact `var.feature_name` match:

- LOX: ENSG00000113083
- LOXL1: ENSG00000129038
- LOXL2: ENSG00000134013
- LOXL3: ENSG00000115318
- LOXL4: ENSG00000138131

## Donor-Aware Summary Readiness

Donor-aware pseudobulk input preparation appears feasible. The Stage 3A input summary produced 156 donor/sample/development/sex/subtype cell-count groups and 780 target-gene detection preview rows.

The detection preview uses `X > 0` only. No mean-expression interpretation, statistical test, cross-dataset comparison, or mouse-human comparison was performed.

## Stage 3B Recommendation

- Stage 3B expression summary should proceed: yes, after confirming matrix semantics and deciding whether `X` or `raw.X` is the appropriate matrix for the planned summary.
- Main blockers: matrix/layer semantics remain partial; `celltypes` and `cell_type` should be documented as fine versus broad labels before reporting subtype summaries.

This Stage 3A audit establishes local parsing feasibility and donor-aware summary readiness for the Park human TEC dataset. It does not establish human conservation, validation, mechanism, protein-level evidence, or therapeutic relevance.
