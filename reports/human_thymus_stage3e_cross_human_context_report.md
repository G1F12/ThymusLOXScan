# Stage 3E Cross-Human LOX-Family Context Report

Search date: 2026-07-06

## Purpose

Stage 3E builds a cautious cross-human context matrix from prior donor/sample-aware or sample-aware Stage 3 outputs. It uses prior summary tables only and does not run new cell-level or cross-dataset statistical tests.

## Datasets Included

- Park TEC: donor/sample-aware human developmental TEC context.
- Yayon TEC: donor/sample-aware human TEC developmental or early-life context.
- GSE147520 epithelial: sample-aware human epithelial context; donor and sex fields unavailable in the inspected H5AD.

## Dataset-Level Overview

Park and Yayon have donor, sample, development-stage, sex, fine annotation, and broad annotation fields. GSE147520 has sample and epithelial-label fields, with age/development information embedded in sample labels.

## Gene-Level LOX-Family Coverage

All five target genes have detection-summary context in the prior outputs. Park and Yayon use X for detection. GSE147520 uses raw.X because compact X omits the targets.

## Donor/Sample-Aware Status

Park and Yayon are donor/sample-aware. GSE147520 is sample-aware only because donor fields were unavailable in the inspected H5AD.

## Matrix-Semantics Status

All three datasets remain partial for matrix value semantics. Mean-value synthesis is therefore not generated.

## Age/Development Coverage

Park and Yayon mainly support developmental or early-life TEC context. GSE147520 includes fetal, postnatal, and one adult-labeled sample group, but adult coverage remains limited.

## What the Three Datasets Support

Together, these outputs support exploratory donor/sample-aware or sample-aware human TEC/epithelial detection-context summaries for LOX-family genes.

## What They Do Not Support

They do not establish aged-adult human thymus behavior, do not support protein or functional interpretation, and do not support cross-species conclusions.

## GSE231906 Next

GSE231906 should be pursued next only as a controlled large-data pilot if aged-adult human relevance is the priority. Its archive size and join structure remain the main blockers.

## Manuscript Decision

Manuscript update is not recommended now. The current result is useful internal context but should remain outside the manuscript until aged-adult or orthogonal follow-up is available.

Across Park TEC, Yayon TEC, and GSE147520 epithelial outputs, the human analyses provide exploratory transcript-level TEC/epithelial context for LOX-family genes. They do not establish human conservation of the mouse mTEC1 Loxl2 candidate signal, do not validate the mouse result, and do not provide mechanism, protein-level evidence, functional evidence, LOX activity evidence, or therapeutic relevance.

Figures created: results\figures\human_thymus_stage3e_cross_context\stage3e_dataset_coverage_overview.png, results\figures\human_thymus_stage3e_cross_context\stage3e_lox_gene_detection_coverage.png, results\figures\human_thymus_stage3e_cross_context\stage3e_dataset_limitations_overview.png