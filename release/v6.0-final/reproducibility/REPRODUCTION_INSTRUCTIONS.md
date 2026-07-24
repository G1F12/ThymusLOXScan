# Reproducibility notes

This package documents frozen outputs rather than a newly executed analysis. Public inputs can be obtained from GEO (GSE240016, GSE223049, and GSE231906) and BioStudies (E-MTAB-8560). See `../governance/datasets_v6.tsv` and `input_files_v6.tsv` for accessions and input provenance.

Use `environment.yml` and `requirements.txt` as the supplied environment specifications. Frozen compact tables can be checked directly against the final manuscript, figures, and checksums index.

The controlled consolidation command is:

`python scripts/run_phase9_final_computational_freeze.py`

Two isolated `FULL_LOCAL` attempts reproduced the available Python core but remained partial because specified external inputs and the R/Bioconductor environment were unavailable. `RELEASE_DERIVED` checks completed without numerical invariant failures. The release therefore does not claim complete clean-room reproduction of every historical workflow.
