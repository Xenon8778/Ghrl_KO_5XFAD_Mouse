library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(AUCell)
library(patchwork)
library(ComplexHeatmap)
library(dittoSeq)
library(scales)
library(openxlsx)
library(RColorBrewer)

# Color maps — adapt as needed
colmap_genotype <- c('5XFAD' = 'steelblue1', '5XFAD; Ghrl -/-' = 'goldenrod1')

#--------------
# Load & subset Microglia
#--------------
so <- readRDS('data/Merged_CellBender_Ghrelin_annotated.rds')
so$ghrelin  <- factor(so$ghrelin,  levels = c('WT', 'Ghrl KO'))
so$genotype <- recode_factor(so$ghrelin,
                             'WT'      = '5XFAD',
                             'Ghrl KO' = '5XFAD; Ghrl -/-')

mg <- subset(so, subset = celltype == 'Microglia')
mg <- NormalizeData(mg)
DefaultAssay(mg) <- 'RNA'

#==============================================================
# 1. AUCell Scoring
#==============================================================

# Gene signatures
signatures <- list(
  DAM         = c('Cstb', 'Fth1', 'Trem2', 'Irf8', 'Timp2', 'Tyrobp',
                  'Ccl6', 'Ccl3', 'Lyz2', 'Axl', 'Cst7', 'B2m', 'Cd9',
                  'Lpl', 'Csf1', 'Cd63', 'Clec7a', 'Igf1', 'Ctsd', 'Lgals3',
                  'Itgax'),
  Homeostatic = c('Jun', 'Hexb', 'Ms4a6d', 'Bin1', 'Sall1', 'Selplg', 'Siglech',
                  'Mef2a', 'Serinc3', 'Tgfbr1',
                  'P2ry12', 'Cd33', 'Tmem119', 'Csf1r', 'Cx3cr1')
)

# Build expression matrix and rank
expr_mat <- GetAssayData(mg, assay = 'RNA', layer = 'counts')
cell_rankings <- AUCell_buildRankings(expr_mat, plotStats = FALSE, verbose = FALSE)

# Calculate AUC
auc_scores <- AUCell_calcAUC(signatures, cell_rankings, verbose = FALSE)

# Add scores to metadata
mg$AUC_DAM         <- as.numeric(getAUC(auc_scores)['DAM', ])
mg$AUC_Homeostatic <- as.numeric(getAUC(auc_scores)['Homeostatic', ])

# ── Violin plots ────────────────────────────────────────────
comparisons <- list(c('5XFAD', '5XFAD; Ghrl -/-'))

p_dam <- ggplot(mg@meta.data, aes(x = genotype, y = AUC_DAM, fill = genotype)) +
  geom_violin(trim = FALSE, scale = 'width') +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = 'white', alpha = 0.6) +
  geom_signif(comparisons = comparisons,
              map_signif_level = TRUE, margin_top = 0.2, vjust = 2) +
  scale_fill_manual(values = colmap_genotype) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = 'DAM Signature', x = NULL, y = 'AUCell Score')

p_hom <- ggplot(mg@meta.data, aes(x = genotype, y = AUC_Homeostatic, fill = genotype)) +
  geom_violin(trim = FALSE, scale = 'width') +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = 'white', alpha = 0.6) +
  geom_signif(comparisons = comparisons,
              map_signif_level = TRUE, margin_top = 0.2, vjust = 2) +
  scale_fill_manual(values = colmap_genotype) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = 'Homeostatic Signature', x = NULL, y = 'AUCell Score')

p_aucell <- p_dam + p_hom +
  plot_annotation(title = 'Microglia AUCell Scores by Genotype')
p_aucell

png('plots/Microglia_AUCell_Violin.png', width = 7, height = 5, units = 'in', res = 300)
p_aucell
dev.off()

#==============================================================
# 1b. AddModuleScore Scoring
#==============================================================
mg <- AddModuleScore(mg,
                     features = list(signatures$DAM),
                     name     = 'ModScore_DAM')
mg <- AddModuleScore(mg,
                     features = list(signatures$Homeostatic),
                     name     = 'ModScore_Homeostatic')

# Seurat appends '1' to the name
colnames(mg@meta.data)[colnames(mg@meta.data) == 'ModScore_DAM1']         <- 'ModScore_DAM'
colnames(mg@meta.data)[colnames(mg@meta.data) == 'ModScore_Homeostatic1'] <- 'ModScore_Homeostatic'

# ── Violin plots ────────────────────────────────────────────
p_mod_dam <- ggplot(mg@meta.data, aes(x = genotype, y = ModScore_DAM, fill = genotype)) +
  geom_violin(trim = TRUE, scale = 'width', alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = 'white', alpha = 0.6) +
  geom_signif(comparisons = comparisons,
              map_signif_level = TRUE, margin_top = 0.05, vjust = 2) +
  scale_fill_manual(values = colmap_genotype) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = 'DAM Signature', x = NULL, y = 'Module Score')

p_mod_hom <- ggplot(mg@meta.data, aes(x = genotype, y = ModScore_Homeostatic, fill = genotype)) +
  geom_violin(trim = TRUE, scale = 'width', alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = 'white', alpha = 0.6) +
  geom_signif(comparisons = comparisons,
              map_signif_level = TRUE, margin_top = 0.05, vjust = 2) +
  scale_fill_manual(values = colmap_genotype) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = 'Homeostatic Signature', x = NULL, y = 'Module Score')

p_modscore <- p_mod_dam + p_mod_hom +
  plot_annotation(title = 'Microglia Module Scores by Genotype')
p_modscore

png('plots/Microglia_ModuleScore_Violin.png', width = 7, height = 5, units = 'in', res = 300)
p_modscore
dev.off()

## DAM Heatmap
avg_expr <- AverageExpression(
  mg, layer = 'data',
  assays        = 'RNA',
  group.by      = c('genotype','sample'),
  return.seurat = FALSE
)$RNA[unique(signatures$DAM),]

# Parse column metadata
col_meta <- data.frame(
  col_id   = colnames(avg_expr),
  genotype = substring(sapply(strsplit(colnames(avg_expr), '_'), `[`, 1),2),   # first field
  sample   = sapply(strsplit(colnames(avg_expr), '_'), `[`, 2),
  stringsAsFactors = FALSE
)

# Sort columns: genotype → sample
col_meta <- col_meta[order(col_meta$genotype,
                           col_meta$sample), ]
avg_expr <- avg_expr[, col_meta$col_id]

# Scale rows for visibility
avg_scaled <- t(scale(t(avg_expr)))
avg_scaled[avg_scaled >  3] <-  3
avg_scaled[avg_scaled < -3] <- -3

# ── Colour palettes ──────────────────────────────────────────
geno_cols    <- colmap_genotype   # reuse project palette

# ── Column annotation: Genotype (top) + Sample (bottom) ─────
ha_col <- HeatmapAnnotation(
  Genotype = col_meta$genotype,
  col      = list(Genotype = geno_cols),
  annotation_name_side = 'right'
)

ht <- Heatmap(
  avg_scaled,
  name                  = 'Scaled\nExpr',
  col                   = circlize::colorRamp2(c(-3, 0, 3),
                                               c('dodgerblue3', 'white', 'firebrick3')),
  top_annotation        = ha_col,
  cluster_rows          = FALSE,          # rows ordered by cluster
  cluster_columns       = FALSE,          # columns ordered by genotype→sample
  column_split          = col_meta$cluster,
  row_title_gp          = gpar(fontsize = 8, fontface = 'bold'),
  column_title_gp       = gpar(fontsize = 8, fontface = 'bold'),
  row_names_gp          = gpar(fontsize = 10),
  column_names_rot      = 90,
  show_column_names     = TRUE,
  use_raster            = TRUE
)
draw(ht)

png('plots/Microglia_DAM_MarkerHeatmap.png',
    width = 4, height = 6, units = 'in', res = 300)
draw(ht)
dev.off()

## Homeostatic Heatmap
avg_expr <- AverageExpression(
  mg, layer = 'data',
  assays        = 'RNA',
  group.by      = c('genotype','sample'),
  return.seurat = FALSE
)$RNA[unique(signatures$Homeostatic),]

# Parse column metadata
col_meta <- data.frame(
  col_id   = colnames(avg_expr),
  genotype = substring(sapply(strsplit(colnames(avg_expr), '_'), `[`, 1),2),   # first field
  sample   = sapply(strsplit(colnames(avg_expr), '_'), `[`, 2),
  stringsAsFactors = FALSE
)

# Sort columns: genotype → sample
col_meta <- col_meta[order(col_meta$genotype,
                           col_meta$sample), ]
avg_expr <- avg_expr[, col_meta$col_id]

# Scale rows for visibility
avg_scaled <- t(scale(t(avg_expr)))
avg_scaled[avg_scaled >  3] <-  3
avg_scaled[avg_scaled < -3] <- -3

# ── Colour palettes ──────────────────────────────────────────
geno_cols    <- colmap_genotype   # reuse project palette

# ── Column annotation: Genotype (top) + Sample (bottom) ─────
ha_col <- HeatmapAnnotation(
  Genotype = col_meta$genotype,
  col      = list(Genotype = geno_cols),
  annotation_name_side = 'right'
)

ht <- Heatmap(
  avg_scaled,
  name                  = 'Scaled\nExpr',
  col                   = circlize::colorRamp2(c(-3, 0, 3),
                                               c('dodgerblue3', 'white', 'firebrick3')),
  top_annotation        = ha_col,
  cluster_rows          = FALSE,          # rows ordered by cluster
  cluster_columns       = FALSE,          # columns ordered by genotype→sample
  column_split          = col_meta$cluster,
  row_title_gp          = gpar(fontsize = 8, fontface = 'bold'),
  column_title_gp       = gpar(fontsize = 8, fontface = 'bold'),
  row_names_gp          = gpar(fontsize = 10),
  column_names_rot      = 90,
  show_column_names     = TRUE,
  use_raster            = TRUE
)
draw(ht)

png('plots/Microglia_Homeostatic_MarkerHeatmap.png',
    width = 4, height = 5, units = 'in', res = 300)
draw(ht)
dev.off()


#==============================================================
# 2. Microglia Subclustering + Marker Heatmap
#==============================================================

set.seed(42)
mg <- FindVariableFeatures(mg, nfeatures = 2000)
mg <- ScaleData(mg)
mg <- RunPCA(mg, npcs = 30)
mg <- RunUMAP(mg, dims = 1:20)
mg <- FindNeighbors(mg, dims = 1:20)
mg <- FindClusters(mg, resolution = 0.2, algorithm = 'leiden',
                   cluster.name = 'mg_subcluster', random.seed = 42)

# ── Save annotated object ────────────────────────────────────
saveRDS(mg, 'data/Microglia_subclustered.rds')
message('Microglia analysis complete.')

#==============================================================
# Panel — UMAP coloured by subcluster
#==============================================================
mg <- readRDS('data/Microglia_subclustered.rds')

p_subclust <- DimPlot(mg, group.by = 'mg_subcluster', label = TRUE, repel = TRUE) +
  ggtitle('Microglia Subclusters') + NoLegend()+
  coord_fixed()
p_subclust
png('plots/Microglia_Subcluster_UMAP.png',
    width = 5, height = 4, units = 'in', res = 300)
print(p_subclust)
dev.off()

# ── Marker genes ────────────────────────────────────────────
Idents(mg) <- 'mg_subcluster'
mg_markers <- FindAllMarkers(mg, assay = 'RNA', only.pos = TRUE,
                             logfc.threshold = 1, min.pct = 0.1) %>%
  filter(p_val_adj < 0.05)

top10 <- mg_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(mg_markers, 'results/Microglia_subcluster_markers.csv', row.names = FALSE)

# ── Average expression per genotype × sample × cluster ──────
# Key: group by genotype__sample__cluster so columns are nested
# genotype > sample, and each cluster gets its own column block
mg$genotype_sample_cluster <- paste0(mg$genotype, '::', mg$sample, '::', mg$mg_subcluster)

avg_expr <- AverageExpression(
  mg, layer = 'data',
  assays        = 'RNA',
  group.by      = c('genotype','sample','mg_subcluster'),
  return.seurat = FALSE
)$RNA[unique(top10$gene),]

# Parse column metadata
col_meta <- data.frame(
  col_id   = colnames(avg_expr),
  genotype = substring(sapply(strsplit(colnames(avg_expr), '_'), `[`, 1),2),   # first field
  sample   = sapply(strsplit(colnames(avg_expr), '_'), `[`, 2),
  cluster  = sapply(strsplit(colnames(avg_expr), '_'), `[`, 3),
  stringsAsFactors = FALSE
)

# Sort columns: cluster → genotype → sample
col_meta <- col_meta[order(as.integer(col_meta$cluster),
                           col_meta$genotype,
                           col_meta$sample), ]
avg_expr <- avg_expr[, col_meta$col_id]

# Row order: genes grouped by their source cluster
gene_order <- unique(top10$gene[order(as.integer(as.character(top10$cluster)))])
gene_order <- intersect(gene_order, rownames(avg_expr))
avg_expr   <- avg_expr[gene_order, ]

# Scale rows for visibility
avg_scaled <- t(scale(t(avg_expr)))
avg_scaled[avg_scaled >  3] <-  3
avg_scaled[avg_scaled < -3] <- -3

# ── Colour palettes ──────────────────────────────────────────
n_clusters   <- length(unique(col_meta$cluster))
clust_cols   <- setNames(hue_pal()(n_clusters), sort(unique(col_meta$cluster)))
geno_cols    <- colmap_genotype   # reuse project palette

# ── Column annotation: Genotype (top) + Sample (bottom) ─────
ha_col <- HeatmapAnnotation(
  Genotype = col_meta$genotype,
  # Sample   = col_meta$sample,
  col      = list(Genotype = geno_cols),
  annotation_name_side = 'left'
)

# ── Row annotation: which cluster each gene belongs to ───────
row_cluster <- top10$cluster[match(gene_order, top10$gene)]

ra_row <- rowAnnotation(
  Cluster = as.character(row_cluster),
  col     = list(Cluster = clust_cols),
  show_annotation_name = FALSE
)

ht <- Heatmap(
  avg_scaled,
  name                  = 'Scaled\nExpr',
  col                   = circlize::colorRamp2(c(-3, 0, 3),
                                               c('dodgerblue3', 'white', 'firebrick3')),
  top_annotation        = ha_col,
  left_annotation       = ra_row,
  cluster_rows          = FALSE,          # rows ordered by cluster
  cluster_columns       = FALSE,          # columns ordered by cluster→genotype→sample
  column_split          = col_meta$cluster,
  row_split             = as.character(row_cluster),
  row_title_gp          = gpar(fontsize = 8, fontface = 'bold'),
  column_title_gp       = gpar(fontsize = 8, fontface = 'bold'),
  row_names_gp          = gpar(fontsize = 8),
  column_names_gp       = gpar(fontsize = 6),
  column_names_rot      = 45,
  show_column_names     = FALSE,
  use_raster            = TRUE
)
draw(ht)

png('plots/Microglia_Subcluster_MarkerHeatmap.png',
    width = 8, height = 8, units = 'in', res = 300)
draw(ht)
dev.off()

#==============================================================
# 3. Proportions of genotype in each subcluster
#==============================================================

# Count table: subcluster × genotype × sample (for error bars)
prop_df <- as.data.frame(table(
  subcluster = mg$mg_subcluster,
  genotype   = mg$genotype,
  sample     = mg$sample
))

prop_df <- prop_df %>%
  group_by(genotype, sample) %>%
  mutate(pct = Freq / sum(Freq) * 100) %>%
  ungroup()

# ── Stacked bar (overall proportions) ───────────────────────
p_stack <- dittoBarPlot(mg, var = 'genotype', group.by = 'mg_subcluster',
                        retain.factor.levels = FALSE,
                        color.panel = unname(colmap_genotype)) +
  geom_col(color = 'black') +
  theme_classic() +
  labs(title = 'Proportion per Microglia Subcluster',
       x = 'Subcluster', y = 'Proportion', fill = 'Genotype') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_stack

png('plots/Microglia_GenotypeProportions_Stacked.png',
    width = 4, height = 3, units = 'in', res = 300)
p_stack
dev.off()

p_h <- dittoBarPlot(mg, var = 'mg_subcluster', group.by = 'genotype', retain.factor.levels = T,
                    color.panel = hue_pal()(6))+
  geom_col(color = 'black') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(title=NULL, x=NULL)
p_h
png('plots/Microglia_Composition_Bar.png',
    width = 3, height = 4, units = 'in', res = 300)
print(p_h)
dev.off()

# ── Save proportion table ────────────────────────────────────
wb <- createWorkbook()
addWorksheet(wb, 'SubclusterProportions', gridLines = FALSE)
writeDataTable(wb, 'SubclusterProportions', prop_df)
saveWorkbook(wb, 'results/Microglia_SubclusterProportions.xlsx', overwrite = TRUE)


# Bar plots by replicate
countDF <- table(mg$mg_subcluster, mg$genotype, mg$sample)
countDF <- reshape2::melt(countDF)
colnames(countDF) <- c('mg_subcluster','genotype','sample','count')
countDF

percDF <- countDF %>%
  group_by(genotype, sample) %>%
  mutate(Percentage = (count / sum(count)) * 100) %>%
  ungroup()

p1 = percDF %>% ggplot(aes(x=genotype, y=Percentage, fill=genotype))+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)+
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7)+
  geom_jitter(aes(shape = genotype), width = 0.2)+
  geom_signif(comparisons = list(c('5XFAD', '5XFAD; Ghrl -/-')),
              test = 't.test',
              map_signif_level = T,margin_top = 0.8, vjust = 1.5)+
  facet_wrap(~mg_subcluster, ncol=6, scales='free')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = colmap_genotype)+
  xlab(NULL)
p1

png('plots/Microglia_Composition_Bar_Signif.png', width = 10, height = 3, units = 'in', res = 300)
p1
dev.off()

#==============================================================
# 4. Pathway AUCell + Heatmaps + DE Overlap Dot Plot
#==============================================================

# ── 4a. Pathway gene signatures ─────────────────────────────
# Mouse gene symbols (sentence-case)
pathway_sigs <- list(
  IGF1_PI3K_AKT_mTOR = c(
    'Igf1', 'Igf1r', 'Irs1', 'Irs2',
    'Pik3ca', 'Pik3cb', 'Pik3cd', 'Pik3r1', 'Pik3r2',
    'Pten',
    'Akt1', 'Akt2', 'Akt3',
    'Mtor', 'Rptor', 'Rictor',
    'S6k1', 'Rps6kb1', 'Eif4ebp1',
    'Foxo1', 'Foxo3',
    'Gsk3b', 'Tsc1', 'Tsc2'
  ),
  TGFb_Signaling = c(
    'Tgfb1', 'Tgfb2', 'Tgfb3',
    'Tgfbr1', 'Tgfbr2', 'Tgfbr3',
    'Smad2', 'Smad3', 'Smad4', 'Smad6', 'Smad7',
    'Ltbp1',
    'Thbs1',
    'Acvr1', 'Acvr2a', 'Acvr2b',
    'Bmpr1a', 'Bmpr1b', 'Bmpr2',
    'Bmp2', 'Bmp4', 'Bmp7',
    'Inhba', 'Inhbb',
    'Ndr1', 'Ndr2',
    'Skil', 'Smurf1', 'Smurf2'
  )
)

# ── 4b. AUCell scoring ───────────────────────────────────────
# Reuse existing cell_rankings if already in environment; rebuild if not
if (!exists('cell_rankings')) {
  expr_mat      <- GetAssayData(mg, assay = 'RNA', layer = 'counts')
  cell_rankings <- AUCell_buildRankings(expr_mat, plotStats = FALSE, verbose = FALSE)
}

pathway_auc <- AUCell_calcAUC(pathway_sigs, cell_rankings, verbose = FALSE)

mg$AUC_IGF1_PI3K  <- as.numeric(getAUC(pathway_auc)['IGF1_PI3K_AKT_mTOR', ])
mg$AUC_TGFb       <- as.numeric(getAUC(pathway_auc)['TGFb_Signaling', ])

# ── Violin plots ─────────────────────────────────────────────
comparisons <- list(c('5XFAD', '5XFAD; Ghrl -/-'))
p_igf <- ggplot(mg@meta.data, aes(x = genotype, y = AUC_IGF1_PI3K, fill = genotype)) +
  geom_violin(trim = FALSE, scale = 'width') +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = 'white', alpha = 0.6) +
  geom_signif(comparisons = comparisons,
              map_signif_level = TRUE, margin_top = 0.2, vjust = 2) +
  scale_fill_manual(values = colmap_genotype) +
  theme_classic() +
  facet_wrap(~mg_subcluster)+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = 'IGF1 / PI3K-AKT-mTOR', x = NULL, y = 'AUCell Score')
p_igf

p_tgf <- ggplot(mg@meta.data, aes(x = genotype, y = AUC_TGFb, fill = genotype)) +
  geom_violin(trim = FALSE, scale = 'width') +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = 'white', alpha = 0.6) +
  geom_signif(comparisons = comparisons,
              map_signif_level = TRUE, margin_top = 0.2, vjust = 2) +
  scale_fill_manual(values = colmap_genotype) +
  theme_classic() +
  facet_wrap(~mg_subcluster)+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = 'TGF-β Receptor Signaling', x = NULL, y = 'AUCell Score')
p_tgf

p_pathway_aucell <- p_igf + p_tgf +
  plot_annotation(title = 'Microglia Pathway AUCell Scores by Genotype')
p_pathway_aucell

png('plots/Microglia_IGF_TGFB_AUCell_Violin.png', width = 10, height = 6, units = 'in', res = 300)
p_pathway_aucell
dev.off()

# ── Helper: build per-sample heatmap for a gene set ──────────
make_pathway_heatmap <- function(mg_obj, genes, title_label, colmap,
                                 group.by   = 'genotype',
                                 ident.1    = '5XFAD; Ghrl -/-',
                                 ident.2    = '5XFAD',
                                 fc_thresh  = 1,
                                 pval_thresh = 0.05) {

  # Keep only genes present in the assay
  genes <- intersect(genes, rownames(mg_obj))

  # ── DE on this gene set ──────────────────────────────────────
  de_res <- tryCatch(
    FindMarkers(mg_obj, features = genes, group.by = group.by,
      ident.1 = ident.1, ident.2 = ident.2,
      logfc.threshold = 0, min.pct = 0.05,test.use = 'wilcox'
    ) %>% mutate(gene = rownames(.)),
    error = function(e) {
      message('DE skipped: ', conditionMessage(e))
      data.frame(gene = character(), p_val_adj = numeric(),
                 avg_log2FC = numeric())
    }
  )

  sig_genes <- de_res %>%
    filter(p_val < pval_thresh, abs(avg_log2FC) > fc_thresh) %>%
    pull(gene)
  print(sig_genes)

  # ── Average expression per genotype × sample ────────────────
  avg_expr <- AverageExpression(
    mg_obj, layer = 'data',
    assays        = 'RNA',
    group.by      = c('genotype', 'sample'),
    return.seurat = FALSE
  )$RNA[genes, , drop = FALSE]

  # Parse column metadata
  split_cols <- strsplit(colnames(avg_expr), '_')
  col_meta <- data.frame(
    col_id   = colnames(avg_expr),
    genotype = substring(sapply(split_cols, `[`, 1), 2),
    sample   = sapply(split_cols, `[`, 2),
    stringsAsFactors = FALSE
  )
  col_meta <- col_meta[order(col_meta$genotype, col_meta$sample), ]
  avg_expr <- avg_expr[genes[genes %in% rownames(avg_expr)],
                       col_meta$col_id, drop = FALSE]

  # Scale rows
  avg_scaled <- t(scale(t(avg_expr)))
  avg_scaled[avg_scaled >  3] <-  3
  avg_scaled[avg_scaled < -3] <- -3

  # ── Row label styling: bold + asterisk for sig DE genes ──────
  row_genes  <- rownames(avg_scaled)
  row_labels <- ifelse(row_genes %in% sig_genes,
                       paste0(row_genes, ' *'),
                       row_genes)
  row_fontface <- ifelse(row_genes %in% sig_genes, 'bold.italic', 'plain')

  # ── Column annotation ────────────────────────────────────────
  ha_col <- HeatmapAnnotation(
    Genotype = col_meta$genotype,
    col      = list(Genotype = colmap),
    annotation_name_side = 'right'
  )

  Heatmap(
    avg_scaled,
    name               = 'Scaled\nExpr',
    col                = circlize::colorRamp2(c(-3, 0, 3),
                                              c('dodgerblue3', 'white', 'firebrick3')),
    top_annotation     = ha_col,
    cluster_rows       = TRUE,
    cluster_columns    = FALSE,
    column_split       = col_meta$genotype,
    row_labels         = row_labels,
    row_names_gp       = gpar(fontsize = 10, fontface = row_fontface),
    column_names_rot   = 90,
    show_column_names  = TRUE,
    column_title       = title_label,
    column_title_gp    = gpar(fontsize = 12, fontface = 'bold'),
    use_raster         = TRUE
  )
}

# ── 4c. Per-gene heatmaps across all samples and genotypes ───
ht_igf <- make_pathway_heatmap(
  subset(mg, subset = mg_subcluster == "1"),
  pathway_sigs$IGF1_PI3K_AKT_mTOR,
  'Sub Cluster 1 - IGF1 / IRS1-PI3K-AKT-mTOR', colmap_genotype,
  fc_thresh  = 0, pval_thresh = 0.05
)
draw(ht_igf)
png('plots/Microglia_SubCluster1_IGF1_PI3K_Heatmap.png',
    width = 6, height = 8, units = 'in', res = 300)
draw(ht_igf)
dev.off()

ht_tgf <- make_pathway_heatmap(
  mg, pathway_sigs$TGFb_Signaling,
  'Sub Cluster 1 - TGF-β Receptor Signaling', colmap_genotype,
  fc_thresh  = 0, pval_thresh = 0.05
)
png('plots/Microglia_SubCluster1_TGFb_Heatmap.png',
    width = 6, height = 9, units = 'in', res = 300)
draw(ht_tgf)
dev.off()

message('Section 4 complete.')


