# GSE240016 mTEC1 sample metadata and batch audit

## Scope

This audit inspected the local processed AnnData file `data/raw/GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad` and local data directories. No local GEO/SRA metadata table for GSE240016 was found.

## Available local metadata

The AnnData object contains `sample`, `stage`, `day`, cell-level QC metrics, and cell annotations. All cells have `day=d0`. It does not expose mouse ID, sequencing run, accession, library preparation date, processing batch, sex, or litter.

Four biological sample labels have at least 10 annotated `13:mTEC1` cells and were used for the mTEC1 Loxl2 sample-level analyses:

- `mo02_CD45neg1_d0`, 02mo, 300 mTEC1 cells
- `mo02_CD45neg2_d0`, 02mo, 279 mTEC1 cells
- `mo18_CD45neg1_d0`, 18mo, 289 mTEC1 cells
- `mo18_CD45neg2_d0`, 18mo, 764 mTEC1 cells

Two aged non-CD45neg-enriched sample labels have only 3 and 2 annotated mTEC1 cells and are excluded by the same minimum-cell threshold used in prior audits.

## Depth and QC

Aged mTEC1 samples have lower median total counts and detected genes than young mTEC1 samples:

- Young mTEC1 median total counts: 8038.5 and 9032.0
- Aged mTEC1 median total counts: 5641.0 and 5457.0
- Young mTEC1 median detected genes: 2819.5 and 3006.0
- Aged mTEC1 median detected genes: 2317.0 and 2175.5

This supports retaining depth/dropout as a residual concern.

## Batch confounding conclusion

Available metadata do not allow batch confounding to be ruled out. Because age and biological sample are highly coupled in the n=2 versus n=2 design, residual sample/batch confounding remains possible.

The local object cannot answer whether age perfectly coincides with sequencing run, library preparation date, or processing batch because those fields are unavailable locally. A fetch step for GEO/SRA sample and run metadata is needed before any stronger batch statement.

## Suggested metadata fetch

Fetch and archive a small derived metadata table, not raw sequence files:

```bash
# Example approach; exact accessions should be confirmed from the GSE240016 GEO record.
python -m pip install pysradb
pysradb metadata GSE240016 --detailed > data/external/GSE240016_metadata.tsv
```

After fetching, join the GEO/SRA metadata to the eight local `sample` labels and test whether age is confounded with library, run, or processing fields.

## Output

- `results/tables/gse240016_sample_metadata_batch_audit.tsv`
