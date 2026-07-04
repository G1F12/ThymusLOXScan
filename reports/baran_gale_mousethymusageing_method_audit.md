# Baran-Gale / MouseThymusAgeing method audit

## Classification

Baran-Gale methodological classification: **provisional directional context: sample identity uncertain**.

True biological mouse/sample aggregation was **not verified** from the cached RDS metadata parsed in the current environment. The current summaries should remain provisional directional TEC/mTEC-like context, not exact GSE240016 mTEC1 confirmation.

## Files used

The current fallback used cached MouseThymusAgeing / ExperimentHub SMART-seq2 RDS files under `data/external/E-MTAB-8560/experimenthub_smartseq_rds/`:

- `rowdata.rds`
- `coldata-day1.rds` through `coldata-day5.rds`
- `counts-processed-day1.rds` through `counts-processed-day5.rds`

The R exporter did not run, and the expected export TSV under `data/external/BaranGale/` was absent.

## Parsed structures

Direct inspection showed `rowdata.rds` and `coldata-day*.rds` parse as Bioconductor DFrame-like objects via Python `rdata`. `counts-processed-day*.rds` parse as data frames via `pyreadr`, with genes as rows and cell IDs as columns.

The parsed metadata fields were `CellID`, `ClusterID`, `Position`, `PlateID`, `Column`, `Row`, `SortType`, `SortDay`, `Age`, `SubType`, and `sizeFactor`.

## Assay or layer

The current fallback used `counts-processed-day*.rds`. These were treated as count-like values and converted to CPM by summing cell columns. The exact assay/layer remains uncertain because the file name indicates processed counts and no R SingleCellExperiment assay metadata was available from the Python fallback. It should not be described as confirmed raw counts.

## Metadata mapping

- Cell identity: `SortType` and `SubType`
- TEC subtype: `SortType` for cTEC/mTEClo/mTEChi/gmTEC; `SubType` for post-Aire and other mTEC-like labels
- Age: `Age`, parsed to 1, 4, 16, 32, and 52 weeks
- Biological mouse/sample: not verified
- Batch/dataset/run: `PlateID` and `SortDay` exist, but no run/accession/batch or mouse ID was verified

## Sample-proxy meaning

The current script assigns `sample_id = sort_day_<SortDay>`. Here, sample-proxy means a SortDay-derived stratum within an age group, used because no donor-level field was exposed in parsed colData. This is not confirmed to be a biological animal.

Sample-proxy counts per age were usually five per age point, but mTEClo at 52 weeks and post-Aire mTEC at 1 and 52 weeks had four sample proxies.

## Replicate handling

Cells are not used as independent replicates in the current by-sample summary tables; cells are first collapsed into SortDay-derived sample-proxy summaries. However, any inferential p-values over those proxies are not valid biological-sample p-values because true animal/sample identity is unresolved.

## Suitability

- Exact confirmation: no
- Subtype-specific confirmation: no
- Broad TEC/mTEC-like directional context: yes, provisional only
- Feasibility/provisional context: yes

The Baran-Gale result should be described as broad directional context only unless an R-based export or metadata source confirms biological mouse/sample aggregation.

## Output

- `results/tables/external_baran_gale_method_audit.tsv`
