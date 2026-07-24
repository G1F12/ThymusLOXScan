# Release Governance Precedence

## Status

This document resolves authority overlap between scientific evidence governance and public release packaging. It does not freeze the repository, publish a release, rerun an analysis, or independently establish a scientific number.

Scientific authority and release-package authority are different concerns. A file can be scientifically authoritative without being a public release asset, and inclusion in a release manifest cannot make a stale or unresolved scientific output authoritative.

## Precedence by question

Use the narrowest applicable machine-readable registry:

1. **Scientific claims:** `manifests/claim_registry_v6.tsv`
2. **Scientific output status:** `manifests/output_authority_v6.tsv`
3. **Dataset identity and permitted scope:** `manifests/datasets_v6.tsv`
4. **Input provenance:** `manifests/input_files_v6.tsv`
5. **Public package inclusion and release checksums:** `RELEASE_MANIFEST.tsv`
6. **Historical repository-state snapshot:** `manifests/repository_version_state_v6.json`, subject to `manifests/repository_version_state_v6_supersession.json`
7. **Human-readable authority index:** `RELEASE_AUTHORITY.md`
8. **Audit and reporting documents:** `PUBLIC_RELEASE_REMEDIATION_REPORT.md` and related reports

## Conflict rules

- `RELEASE_MANIFEST.tsv` controls public packaging only. It cannot upgrade a stale, superseded, mislabeled, unresolved, or missing-path scientific output into an authoritative result.
- `manifests/output_authority_v6.tsv` controls scientific output status. It cannot automatically add a file to the public release.
- `manifests/claim_registry_v6.tsv` takes precedence over overlapping rows in `results/tables/claim_evidence_registry_v5_6.tsv`. The v5.6 registry may remain provisional support only for claims not represented in v6 and may not override v6.
- The v6 claim registry is a partial registry and must not be represented as a complete claim set.
- Historical files remain part of provenance but are not current public authority.
- Human-readable documents summarize and link to machine-readable authority; they do not replace it.
- When registries conflict, apply the registry whose declared scope most narrowly answers the question. Unresolved conflicts remain unresolved rather than being decided by file date, version label, manuscript citation, or release inclusion.

## GSE231906 controlling boundary

GSE231906 is `DETECTION_CONTEXT_ONLY`.

The strict current join values are 387,762 uniquely matched expression barcodes of 590,533, with aggregate exact unit-scoped match rate 0.656631. Their source table is `results/tables/gse231906_barcode_join_audit_v6.tsv`; scientific scope is controlled by the dataset, output, and claim registries.

The older 388,986 of 590,533 and 0.658703 values are `HISTORICAL/SUPERSEDED`. They were produced under permissive substring-based unit filtering and cannot override the strict v6 join authority.

No GSE231906 human validation, replication, conservation, donor/person-level inference, aging model, mechanism, exact mTEC1 replication, or cross-unit mean-expression comparison is authorized. Detection among uniquely joined barcodes within sample/library units is authoritative. Mean expression is descriptive within unit only.

## Repository-state authority

`manifests/repository_version_state_v6.json` is preserved byte-for-byte as a historical generation-time snapshot. Its description of an uncommitted package at parent commit `6264d1f4be2a7a377782f2dc66f7f8bf39c5d421` was true immediately before source package commit `1f2e422a0da6fc303ebf976abba7dbb11b6b57f5` was created.

`manifests/repository_version_state_v6_supersession.json` removes only the original record's current live-state authority. No replacement current-state registry exists because the worktree is dirty. A new live-state registry must be generated only during release freeze and bound to the exact immutable release commit and tag.

## Release status

This hierarchy does not declare the repository ready for release. Clean-room evidence, dependency locking, provenance gaps, incomplete claim/output lineage, historical release cleanup, manifest freeze, and DOCX visual QA remain separate gates.
