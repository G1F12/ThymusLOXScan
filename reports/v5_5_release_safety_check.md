# v5.5 Release Safety Check

## Scope

- Manuscript source: `manuscript/LOX_thymus_aging_public_preprint_v5_5_external_context.md`
- Release PDF: `manuscript/pdf/LOX_thymus_aging_public_preprint_v5_5_external_context.pdf`
- Release HTML: `manuscript/pdf/LOX_thymus_aging_public_preprint_v5_5_external_context.html`
- Release notes: `reports/v5_5_release_notes.md`

## Artifact Checks

- Release PDF exists: yes.
- Release HTML exists: yes.
- Release PDF title is the scientific title: yes.
- Release PDF has no review notice: yes.
- Release HTML has no review notice: yes.
- Release PDF has no browser headers/footers: yes.
- Release PDF has no file URI paths: yes.
- Release PDF has no local Windows drive-letter paths: yes.
- Figures embedded in HTML: yes.
- Inline Figure 1 renders in the PDF: yes.

## Public Pointer Checks

- README now points the current manuscript to v5.5: yes.
- README now points the current release to `v5.5-external-mtec-context`: yes.
- README lists the public PDF path: yes.
- `one_page_summary.md` now points the current manuscript to v5.5: yes.
- `one_page_summary.md` lists current release `v5.5-external-mtec-context`: yes.
- v5.2 is retained as a previous release: yes.

## Scientific Claim Checks

- No validation, confirmation, proof, robust replication, or exact mTEC1 replication claim for E-MTAB-8560: yes.
- No batch/dropout/depth ruled-out claim: yes.
- No protein/function/mechanism/therapy/human-conservation claim: yes.
- E-MTAB-8560 described as mixed/inconclusive: yes.
- GSE240016 mTEC1 `Loxl2` described as candidate-level: yes.
- Cells are not treated as biological replicates: yes.

Phrase-scan hits for `validation`, `mechanism`, `human conservation`, `rejuvenation`, `therapeutic`, `batch ruled out`, `protein-level validation`, and `exact mTEC1 validation` were reviewed. Hits were file names, explicit negations, limitation wording, or safety-report wording.

## Staging Hygiene Checks

- No raw/large data staged: yes.
- No `.h5ad` staged: yes.
- No `data/raw`, `data/processed`, or large `data/external` staged: yes.
- No `private_outreach` files staged: yes.
- No zip files staged: yes.
- No `environment.yml` staged: yes.
- No `requirements.txt` staged: yes.
- No unrelated dirty/untracked files staged: yes.
- No generated-assistant provenance traces in public files: yes.

## Result

Safety passed: yes.

Release/tag created at this checkpoint: no.
