# E-MTAB-8560 Biological Unit Audit

Classification: infeasible in current environment.

Rscript was not available, so the official MouseThymusAgeing R API could not be used to export colData and assay layers. The SDRF metadata can identify candidate fields such as `Characteristics[individual]`, but this does not verify that every analyzed cell in the official processed object maps to a true biological mouse/individual.

No per-mouse inference was performed. SortDay, PlateID, and run fields were not substituted for biological mouse ID.
