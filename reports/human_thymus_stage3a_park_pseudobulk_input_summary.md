# Park Human TEC Pseudobulk Input Summary

Search date: 2026-07-06

The Park TEC H5AD was summarized into donor-aware input tables for a later expression-summary stage.

- Selected donor field: donor_id (selected_preferred)
- Selected sample field: sample_id (selected_preferred)
- Selected age/development field: development_stage (selected_first_preferred_with_alternatives:development_stage_ontology_term_id)
- Selected sex field: sex (selected_preferred)
- Selected cell-type field: celltypes (selected_first_preferred_with_alternatives:cell_type)
- Cell-count groups: 156
- Target genes present: LOX, LOXL1, LOXL2, LOXL3, LOXL4
- Expression preview generated: yes, detection fraction only.
- Matrix/layer semantics: partial; X is present and raw exists, but layer/count semantics were not asserted by this Stage 3A audit; detection fractions use X > 0 only.

No statistical tests were run, and cells were not treated as biological replicates.