# LOX-family detection-rate analysis by sample

## Scope

Detection was computed per biological sample using raw counts in `adata.raw.X`. A cell is counted as detecting a gene when its raw count for that gene is greater than zero. Mean nonzero expression is the mean raw count among detecting cells only.

This is descriptive and sample-level. It is intended to separate changes in the fraction of expressing cells from changes in expression intensity among cells that still express the gene.

## Key finding interpretation

| group | gene | detection 02mo | detection 18mo | nonzero mean 02mo | nonzero mean 18mo | interpretation |
|---|---:|---:|---:|---:|---:|---|
| pooled FB | Lox | 0.301 | 0.109 | 1.784 | 1.379 | both fewer expressing cells and lower expression among expressing cells |
| pooled FB | Loxl1 | 0.684 | 0.486 | 3.596 | 2.414 | both fewer expressing cells and lower expression among expressing cells |
| pooled FB | Loxl2 | 0.282 | 0.151 | 1.312 | 1.236 | fewer expressing cells |
| capsFB | Lox | 0.343 | 0.120 | 1.703 | 1.156 | both fewer expressing cells and lower expression among expressing cells |
| medFB | Loxl1 | 0.171 | 0.219 | 1.591 | 1.387 | unclear for fewer/lower framework; detection increases while nonzero expression does not |
| medFB | Loxl2 | 0.199 | 0.070 | 1.178 | 1.095 | fewer expressing cells |
| mTEC1 | Loxl2 | 0.090 | 0.003 | 1.058 | 1.000 | fewer expressing cells |

## Notes by finding

### Lox in pooled FB

- Detection-rate change, 18mo minus 02mo: -0.192.
- Mean-nonzero-expression change, 18mo minus 02mo: -0.405 raw counts among detecting cells.
- Expected pseudobulk direction: down_in_aged.
- Classification: both fewer expressing cells and lower expression among expressing cells.

### Loxl1 in pooled FB

- Detection-rate change, 18mo minus 02mo: -0.197.
- Mean-nonzero-expression change, 18mo minus 02mo: -1.182 raw counts among detecting cells.
- Expected pseudobulk direction: down_in_aged.
- Classification: both fewer expressing cells and lower expression among expressing cells.

### Loxl2 in pooled FB

- Detection-rate change, 18mo minus 02mo: -0.130.
- Mean-nonzero-expression change, 18mo minus 02mo: -0.075 raw counts among detecting cells.
- Expected pseudobulk direction: down_in_aged.
- Classification: fewer expressing cells.

### Lox in capsFB

- Detection-rate change, 18mo minus 02mo: -0.223.
- Mean-nonzero-expression change, 18mo minus 02mo: -0.547 raw counts among detecting cells.
- Expected pseudobulk direction: down_in_aged.
- Classification: both fewer expressing cells and lower expression among expressing cells.

### Loxl1 in medFB

- Detection-rate change, 18mo minus 02mo: 0.048.
- Mean-nonzero-expression change, 18mo minus 02mo: -0.203 raw counts among detecting cells.
- Expected pseudobulk direction: up_in_aged.
- Classification: unclear for fewer/lower framework; detection increases while nonzero expression does not.

### Loxl2 in medFB

- Detection-rate change, 18mo minus 02mo: -0.129.
- Mean-nonzero-expression change, 18mo minus 02mo: -0.083 raw counts among detecting cells.
- Expected pseudobulk direction: down_in_aged.
- Classification: fewer expressing cells.

### Loxl2 in mTEC1

- Detection-rate change, 18mo minus 02mo: -0.087.
- Mean-nonzero-expression change, 18mo minus 02mo: -0.058 raw counts among detecting cells.
- Expected pseudobulk direction: down_in_aged.
- Classification: fewer expressing cells.

## Caveats

Detection rates are sensitive to sequencing depth, library complexity, and the small number of biological samples. These summaries should be read alongside the pseudobulk age models rather than as standalone inference.