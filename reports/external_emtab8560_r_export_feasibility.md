# E-MTAB-8560 R Export Feasibility

Rscript is not available in the current environment, so the official MouseThymusAgeing R API could not be run.

No Python/RDS fallback was used. Per-mouse E-MTAB-8560 inference was stopped.

Required installation example:

```r
install.packages("BiocManager")
BiocManager::install(c("MouseThymusAgeing", "SingleCellExperiment", "SummarizedExperiment", "scuttle"))
```

No E-MTAB-8560 per-mouse reanalysis was performed.
