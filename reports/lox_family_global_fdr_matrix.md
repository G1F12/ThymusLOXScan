# LOX-family global FDR matrix

## Scope

This report applies a single Benjamini-Hochberg correction across the internal GSE240016 LOX-family pseudobulk matrix only. External datasets are not mixed into this FDR calculation.

Tested rows with raw p-values: 123.

## Global FDR < 0.05

- FB / Lox: log2FC=-1.295, raw p=1.07e-06, global FDR=6.57e-05.
- 3:capsFB / Lox: log2FC=-1.462, raw p=5.38e-07, global FDR=6.57e-05.
- 13:mTEC1 / Loxl2: log2FC=-3.291, raw p=8.52e-05, global FDR=0.0035.
- 5:medFB / Loxl2: log2FC=-0.995, raw p=0.000441, global FDR=0.0136.
- FB / Loxl1: log2FC=-0.616, raw p=0.00117, global FDR=0.0288.
- FB / Loxl2: log2FC=-0.517, raw p=0.00194, global FDR=0.0385.
- 5:medFB / Loxl1: log2FC=0.751, raw p=0.00219, global FDR=0.0385.
- 4:intFB / Loxl1: log2FC=-0.512, raw p=0.0029, global FDR=0.0446.

## Global FDR < 0.10 only

- 4:intFB / Lox: log2FC=-0.626, raw p=0.00436, global FDR=0.0596.

## Not significant after global matrix correction

114 tested rows were not significant at global FDR < 0.10.

## Output

- `results/tables/lox_family_global_fdr_matrix.tsv`