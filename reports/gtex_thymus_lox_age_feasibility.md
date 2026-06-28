# GTEx thymus LOX-family age feasibility

Classification: age analysis not possible from available public/local files

GTEx remains a planned human whole-tissue analysis requiring donor-level expression and age metadata.

## What was checked

| requirement | status | detail |
|---|---|---|
| local_search_paths | checked | D:\ThymusLOXScan\data\external\GTEx; D:\ThymusLOXScan\data\external\gtex |
| candidate_files_found | 0 |  |
| gene_tpm_expression_matrix | missing |  |
| sample_attributes | missing |  |
| subject_phenotype_age_metadata | missing |  |
| classification | age analysis not possible from available public/local files | missing_files |

## Manual download instructions

Place the GTEx v8 gene TPM matrix and annotation files under `data/external/GTEx/` or `data/external/gtex/`.

Expected local files include:

- `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz` or another gene-level TPM/expression matrix.
- `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt`.
- `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt`.

Do not commit the full GTEx matrix or other large third-party data files to this repository.

## Methods-ready note

If donor-level files are added locally, this analysis will subset GTEx whole-thymus samples, summarize LOX, LOXL1, LOXL2, LOXL3, and LOXL4 as log2(TPM+1), and fit descriptive age associations with available covariates. GTEx thymus is whole tissue, not TEC/mTEC subtype-resolved, and cannot validate mouse mTEC1.

## Discussion-ready note

At this stage, GTEx provides no completed age-context result for this repository. Any future GTEx result should be described as broad human thymus whole-tissue directional context only, and only if donor-level expression and age metadata support that description.