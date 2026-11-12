# Ghrl_KO_5XFAD_Mouse
Code for the single-cell RNA sequencing analysis of hippocampus and medial prefrontal cortex from 5XFAD and 5XFAD; Ghrl−/− mice (n = 3 per genotype, 13–14 months). Ambient RNA was removed with CellBender, and samples were integrated with Harmony and clustered with Leiden in Seurat v5. The pipeline covers:

- Cell type annotation and composition analysis
- Per-cell-type differential expression with GO/KEGG enrichment
- Subclustering of microglia, astrocytes, glutamatergic neurons, oligodendrocytes, OPCs, and choroid plexus
- DAM/homeostatic and pathway signature scoring (AUCell)
- Genotype comparison of cell–cell communication (CellChat)
