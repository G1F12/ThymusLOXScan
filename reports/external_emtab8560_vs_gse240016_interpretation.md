# E-MTAB-8560 vs GSE240016 Interpretation

## GSE240016 Summary

The current GSE240016 mTEC1 Loxl2 signal is transcript-level and candidate-level. It remains limited by n=2 young versus n=2 aged samples, limited exact-permutation resolution, lower aged mTEC1 depth/detected genes, and unresolved batch/sample confounding after official GEO/SRA metadata audit.

Matched-gene falsification makes the observed prioritization less generic, but it does not resolve depth or batch concerns.

## E-MTAB-8560 Status

E-MTAB-8560 remains the best identified independent mouse mTEC-focused ageing dataset candidate. Official BioStudies/ArrayExpress SDRF metadata includes candidate individual/mouse, age, cell-type/sort, and run fields.

However, the official MouseThymusAgeing R export could not be run in the current environment because `Rscript` is unavailable. Therefore, the required cell-level processed count export and colData-to-mouse verification could not be completed.

## Interpretation

This does not change confidence in mTEC/TEC Loxl2 as a transcript-level ageing candidate. No E-MTAB-8560 per-mouse LOX-family result was produced, and no p-values or effect directions should be inferred from this environment.

The recommended next step remains to run the official MouseThymusAgeing R export in an environment with R/Bioconductor installed, then perform the planned per-mouse pseudobulk analysis.

No p-values were combined across datasets.

This is not exact GSE240016 mTEC1 replication.
