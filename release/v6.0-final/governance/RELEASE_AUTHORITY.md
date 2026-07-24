# Release Authority and Document Hierarchy

## Status

This file is the human-readable authority index for a future publication release candidate. It does **not** publish or freeze a release, independently establish scientific numbers, or override the current public v5.6 tag. Detailed precedence is defined in `RELEASE_GOVERNANCE_PRECEDENCE.md`.

The machine-readable authorities are:

- scientific claims: `manifests/claim_registry_v6.tsv`;
- scientific output status: `manifests/output_authority_v6.tsv`;
- dataset identity and permitted scope: `manifests/datasets_v6.tsv`;
- input provenance: `manifests/input_files_v6.tsv`;
- public package inclusion and checksums: `RELEASE_MANIFEST.tsv`.

The claim registry is partial. The input registry contains explicit provenance gaps, and the output registry retains unresolved, stale, mislabeled, and missing-path entries as non-authoritative. None of these gaps is silently promoted by this index.

## Single authoritative manuscript

The only authoritative manuscript source is:

`manuscript/LOX_thymus_aging_public_preprint_v5_6_gse231906_human_context.md`

Its PDF and HTML are renderings, not independent manuscripts. They may be attached only if their hashes match the manifest and their text is consistent with the authoritative Markdown source.

## Deterministic document classification

The following ordered rules give every repository document exactly one role. The first matching rule wins; no manual exception is permitted.

| Priority | Path rule | Role | Publication treatment |
|---:|---|---|---|
| 1 | `manuscript/LOX_thymus_aging_public_preprint_v5_6_gse231906_human_context.md` | `AUTHORITATIVE_MANUSCRIPT` | Required source manuscript |
| 2 | `manuscript/pdf/LOX_thymus_aging_public_preprint_v5_6_gse231906_human_context.{pdf,html}`, `manuscript/docx/LOX_thymus_aging_public_preprint_v5_6_gse231906_human_context.docx`, and its `_submission.zip` | `AUTHORITATIVE_RENDERING_OR_ARCHIVE` | Required only after text/hash consistency gate |
| 3 | `RELEASE_GOVERNANCE_PRECEDENCE.md`, `RELEASE_AUTHORITY.md`, `RELEASE_MANIFEST.tsv`, `REPRODUCIBILITY_PACKAGE.md`, `PUBLIC_RELEASE_REMEDIATION_REPORT.md`, `manifests/*_v6.tsv`, `manifests/repository_version_state_v6*.json` | `RELEASE_GOVERNANCE` | Required or candidate governance according to the release manifest |
| 4 | `manuscript/LOX_thymus_aging_public_preprint_v4_2_final_safe.md`, `manuscript/LOX_thymus_aging_public_preprint_v5*.md` except the authoritative source | `SUPERSEDED` | Exclude from release assets and submission uploads |
| 5 | `manuscript/manuscript.*`, `manuscript/*_revised.md`, `manuscript/*_draft.md`, `manuscript/abstract_revised.md`, `manuscript/discussion_revised.md`, `manuscript/methods_detailed.md`, `manuscript/methods.md`, `manuscript/discussion.md`, `manuscript/results_draft.md` | `ARCHIVED` | Exclude from release assets and submission uploads |
| 6 | `manuscript/pdf/*v5_5*`, all `manuscript/docx/*` not matched by rule 2, and `manuscript/tmp/**` | `REMOVE_FROM_RELEASE` | Do not attach, upload, or expose as a release manuscript |
| 7 | `reports/**`, `results_final_v3.md` | `SUPPORTING_AUDIT_NONAUTHORITATIVE` | May remain in repository history; never cite as an authoritative manuscript or attach as a release manuscript |
| 8 | `one_page_summary.md`, `wet_lab_validation_plan.md`, `lab_outreach/**`, `private_outreach/**` | `OUTREACH_OR_PRIVATE_NONRELEASE` | Exclude from publication release and submission uploads |
| 9 | `notebooks/**`, `scripts/**`, `supplementary_tables/**`, `results/**`, `data/**` | `COMPUTATIONAL_OR_EVIDENCE_MATERIAL` | Scientific status is governed by v6 registries; public inclusion is governed by the release manifest |
| 10 | all remaining documentation files | `REPOSITORY_DOCUMENTATION_NONAUTHORITATIVE` | Exclude unless explicitly listed in the manifest |

## Mandatory citation rule

Only the authoritative source path above may be named as the manuscript in README files, release notes, citations, GitHub releases, Zenodo metadata, or bioRxiv submission metadata. Historical and archived files must not be described as current, final, public manuscript, validated evidence, or release assets.

## GSE231906 publication boundary

For any future release asset, the only allowed GSE231906 scope label is `DETECTION_CONTEXT_ONLY`: target-gene presence/detection in matched metadata-expression units with partial matrix semantics. The label excludes human validation, replication, conservation, donor/person-level inference, aging-model inference, mechanistic support, mean-expression claims, and inclusion in cross-dataset direction synthesis.

The authoritative Markdown source, PDF, HTML, README, one-page summary, and v5.6 release notes now use this boundary. Historical GSE231906 scripts, tables, and reports remain non-authoritative technical materials and are excluded from publication claims and release attachments.

The strict join source is `results/tables/gse231906_barcode_join_audit_v6.tsv`: 387,762 uniquely matched expression barcodes of 590,533, aggregate exact unit-scoped match rate 0.656631. The older 388,986 of 590,533 and 0.658703 values are `HISTORICAL/SUPERSEDED`. These numbers are repeated here only as an index; their authority comes from the machine-readable registries and strict source table.

Detection is authoritative. `COUNT_LIKE_PROBABLE` is not proof of raw-count provenance, and mean expression is descriptive within unit only; cross-unit comparison is not authorized.

## Repository-state records

`manifests/repository_version_state_v6.json` is preserved byte-for-byte as `HISTORICAL_GENERATION_TIME_SNAPSHOT`. It describes the pre-commit state immediately before Phase 0 / Phase 6A source commit `1f2e422a0da6fc303ebf976abba7dbb11b6b57f5`.

Its current live-state authority is superseded by `manifests/repository_version_state_v6_supersession.json`. No replacement current-state registry is created in this dirty worktree. Final live-state authority must be generated during release freeze from the exact immutable release commit and tag.

## Unresolved release gates

- Clean-room reproduction evidence remains incomplete.
- Dependency specifications are not fully locked across Python and R.
- Input provenance still has missing hashes or immutable direct URLs.
- The v6 claim graph and output-authority registry contain unresolved or incomplete entries.
- Historical public/release surfaces still require cleanup or explicit exclusion.
- `RELEASE_MANIFEST.tsv` is not frozen to an immutable commit/tag.
- DOCX visual-render QA remains incomplete.

`PUBLIC_RELEASE_REMEDIATION_REPORT.md` is an audit/remediation report, not primary scientific authority. `REPRODUCIBILITY_PACKAGE.md` documents instructions and status; it is not proof that clean-room reproduction has succeeded.
