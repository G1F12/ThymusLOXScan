# V6 final presentation polish report

## Scope and starting state

- Repository: `D:\ThymusLOXScan`
- Branch: `human-relevance-and-wetlab-plan`
- Starting commit: `32ab4f3d3b7981a0460e40185b5345cdcb23b7fa`
- Scope: presentation and release-documentation polish only
- Scientific workflows executed: none
- Raw or large expression matrices opened: none

The existing worktree contained unrelated modified and untracked material before this task. Those files were preserved and excluded from this polish commit.

## Figures

Exactly three reader-facing figures were created in PNG and PDF formats under `results/figures/v6_final/` and copied byte-for-byte into `release/v6.0-final/figures/`.

1. Figure 1 presents the four frozen GSE240016 focal log2 fold-change estimates: capsFB `Lox` -1.456 (3+3, A), medFB `Loxl2` -0.992 (3+3, B+), medFB `Loxl1` +0.752 (3+3, B), and mTEC1 `Loxl2` -3.286 (2+2, C). No confidence intervals were added.
2. Figure 2 presents the four included mTEC1 `Loxl2` sample-label proxies from the frozen sample-summary table committed at `11ff0c7205bebb5f523bc0163bde6c2f59cef0d7`. Detection fraction is displayed directly; log2(CPM+1) is a deterministic display transform of the frozen CPM column. Cells are not presented as replicates and no inferential P value is shown.
3. Figure 3 maps GSE240016 directions to GSE223049 broad context and E-MTAB-8560 related TEC context. Exact-subtype replication is marked not established. GSE231906 appears only in a separated `DETECTION_CONTEXT_ONLY` scope annotation.

## Manuscript presentation changes

- Replaced the publication-facing status line with `Manuscript version: 6.0`.
- Used `Correspondence: Aliaksandr Karatseyeu` because no verified author email was present.
- Separated the repository as code and materials.
- Rewrote the Abstract to 239 words without letter categories while retaining the primary dataset, 22,932-cell count, four focal directions and sample counts, limitations, mixed/contextual external evidence, GSE231906 detection-only role, and hypothesis-generating conclusion.
- Replaced internal phase jargon in the scientific narrative with descriptive analysis names. Exact phase identifiers remain only in provenance-oriented supplementary mapping.
- Defined A, B+, B, and C once as internal computational robustness categories rather than a universal grading framework.
- Renamed the Table 1 column to `Internal robustness`.
- Removed the full consolidation command from the scientific narrative and retained it only in release/reproducibility/supplementary instructions.
- Added three figure callouts and complete captions identifying source dataset, analytical unit, descriptive or inferential role, and principal limitations.

## Reproducibility wording

Methods, Limitations, Reproducibility and code availability, release README, release notes, and release reproduction instructions now use one interpretation:

> Two isolated `FULL_LOCAL` attempts reproduced the available Python core but remained partial because specified external inputs and the R/Bioconductor environment were unavailable. `RELEASE_DERIVED` checks completed without numerical invariant failures. The release therefore does not claim complete clean-room reproduction of every historical workflow.

This wording does not upgrade clean-room status.

## Rendering and release snapshot

- PDF: 17 A4 pages, numbered in the footer, all pages visually inspected.
- HTML: standalone and self-contained with three embedded figures; no browser headers or local file URLs.
- DOCX: generated as an editable A4 document with three embedded figures and PAGE fields in footers.
- DOCX structural checks passed. A separate DOCX-to-image visual render could not be completed because no Word or LibreOffice-compatible converter was available in the environment; this remains a manual warning.
- Release snapshot: 40 assets, including six figure files.
- Release manifest: one unique row per asset.
- SHA256SUMS: verified for every non-self asset.
- Repository-state record: updated as a presentation-polish snapshot and explicitly left unbound to a final tag or immutable release commit.

## Remaining presentation warning

The clean PDF and HTML were inspected. The DOCX should receive a final visual review in Microsoft Word or LibreOffice before external submission because this environment could not rasterize the DOCX.
