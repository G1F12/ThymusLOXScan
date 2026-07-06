# GSE147520 Human Thymic Epithelial Candidate Context Note

Search date: 2026-07-06

## What GSE147520 Epithelial Can Contribute

The epithelial H5AD provides human thymic epithelial cell-count and subtype context across fetal, postnatal, and one adult-labeled sample group.

## What It Cannot Contribute

The available epithelial H5AD cannot contribute compact-X target-gene summaries because the five target genes are absent from compact X. Detection context can be generated from raw.X, but the file still lacks explicit donor and sex fields in observed metadata.

## Donor-Aware Feasibility

Donor-aware summaries are not fully feasible from this H5AD because a donor field is not available. Sample-aware cell-count context is feasible using `samples`.

## Age and Development Coverage

The `samples` labels include fetal, postnatal, and adult-labeled groups, but adult coverage is limited and does not support direct aged-adult interpretation.

## Workflow Context

GSE147520 complements Park and Yayon at the workflow/context level by showing why matrix-layer inspection matters: compact X omits the targets while raw.X contains them. It should be used cautiously as human epithelial sample-aware context.