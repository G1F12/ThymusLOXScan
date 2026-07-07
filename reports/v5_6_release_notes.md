# v5.6 GSE231906 human aging context pilot

## Summary

- Adds v5.6 public manuscript draft and clean PDF/HTML.
- Adds a controlled GSE231906 human thymus aging pilot.
- Parses 516,352 metadata rows and audits a 3.7 GB processed archive locally without committing large data.
- Identifies 217 archive members and 59 candidate expression units.
- Performs pilot extraction and metadata join for GSM8945751_donor3 with join match rate 0.909477.
- Confirms LOX-family target genes are present: LOX, LOXL1, LOXL2, LOXL3, LOXL4.
- Confirms target compartments are available, including epithelial/Epi, mTEC-like, cTEC-like, post_AIRE_mTEC, immature TEC, mesenchymal/fibroblast-like, and endothelial.
- Generates donor/sample-aware aged-human transcript-level summaries.
- Classifies LOXL2 mTEC-like context as mixed_transcript_level_context.
- Keeps GSE240016 mouse mTEC1 Loxl2 as candidate-level.
- Keeps wet-lab LOXL2 IHC/IF or RNA-ISH as the key next step.

The controlled GSE231906 pilot makes aged-human thymus LOX-family transcript context technically evaluable. LOX-family genes and TEC/mTEC-like compartments were recoverable after metadata-to-expression joining, but the LOXL2 mTEC-like aging context was mixed. This is mixed aged-human transcript-level context, not human validation or conservation of the mouse GSE240016 mTEC1 Loxl2 candidate signal.

## Limitations

- Not human validation.
- Not human conservation.
- Not exact replication.
- Not protein validation.
- Not mechanism.
- Not causality.
- Not functional thymic evidence.
- Not LOX activity evidence.
- Not ECM-crosslinking evidence.
- Not therapeutic evidence.
- Matrix semantics partial.
- Mean-value summaries descriptive only.
- Wet-lab validation still required.
