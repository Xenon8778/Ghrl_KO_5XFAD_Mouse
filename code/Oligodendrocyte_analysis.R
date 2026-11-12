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
library(ggrepel)

colmap_genotype <- c('5XFAD' = 'steelblue1', '5XFAD; Ghrl -/-' = 'goldenrod1')

dir.create('plots',   showWarnings = FALSE)
dir.create('results', showWarnings = FALSE)

#--------------
# Load & subset Oligodendrocytes
#--------------
so <- readRDS('data/Merged_CellBender_Ghrelin_annotated.rds')
so$ghrelin  <- factor(so$ghrelin, levels = c('WT', 'Ghrl KO'))
so$genotype <- recode_factor(so$ghrelin,
                             'WT'      = '5XFAD',
                             'Ghrl KO' = '5XFAD; Ghrl -/-')

oligo <- subset(so, subset = celltype == 'Oligodendrocyte')
oligo <- NormalizeData(oligo)
DefaultAssay(oligo) <- 'RNA'

#==============================================================
# Panel b — Dot plot: key oligodendrocyte marker genes across genotypes
#==============================================================

# Core myelin identity markers (adjust to your biology)
marker_genes_b1 <- c('Mbp','Plp1','Mog','Mag','Cnp')
marker_genes_b1 <- intersect(marker_genes_b1, rownames(oligo))

p_b1 <- DotPlot(oligo,
                features  = marker_genes_b1,
                group.by  = 'genotype', scale.by = 'size', scale.min = 0,
                cols      = c('white', 'firebrick'), scale=F,
                dot.scale = 7) +
  geom_point(aes(size = pct.exp), shape = 21, colour = 'black', stroke = 0.5) +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12)) +
  labs(title = NULL, x = NULL, y = NULL)
p_b1

# Additional maturation / transcription factor markers
marker_genes_b2 <- c('Sox10','Olig1','Olig2','Mobp','Cldn11')
marker_genes_b2 <- intersect(marker_genes_b2, rownames(oligo))

p_b2 <- DotPlot(oligo,
                features  = marker_genes_b2,
                group.by  = 'genotype', scale.by = 'size', scale.min = 0,
                cols      = c('white', 'firebrick'), scale=F,
                dot.scale = 7) +
  geom_point(aes(size = pct.exp), shape = 21, colour = 'black', stroke = 0.5) +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12)) +
  labs(title = NULL, x = NULL, y = NULL)
p_b2

png('plots/Oligodendrocyte_panel_MarkerDotplot.png',
    width = 7, height = 4.5, units = 'in', res = 300)
print(p_b1 + p_b2)
dev.off()

#==============================================================
# Panel d — Dot plot: top recDEGs (upregulated in 5XFAD; Ghrl -/-)
# Panel e — Dot plot: rare recDEGs
# These require DE results — run DE_Analysis.R first
#==============================================================

ident1 <- '5XFAD; Ghrl -/-'
ident2 <- '5XFAD'

# read DE results
WilcoxRNA_sig <- read.xlsx('results/DE_Genotype/GhrlKO_vs_WT_Oligodendrocyte.xlsx', sheet = 'DE_RNA_SIG')
WilcoxRNA_up <- read.xlsx('results/DE_Genotype/GhrlKO_vs_WT_Oligodendrocyte.xlsx', sheet = 'DE_RNA_UP')
WilcoxRNA_down <- read.xlsx('results/DE_Genotype/GhrlKO_vs_WT_Oligodendrocyte.xlsx', sheet = 'DE_RNA_DOWN')

# Panel d — top recDEGs (highest |FC|, expressed in >10% of cells)
top_recDEGs <- WilcoxRNA_sig %>%
  dplyr::filter(pct.1 > 0.1 | pct.2 > 0.1) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n = 10) %>%
  pull(gene)
top_recDEGs <- intersect(top_recDEGs, rownames(oligo))

if (length(top_recDEGs) > 0) {
  p_d <- DotPlot(oligo,
                 features  = top_recDEGs,
                 group.by  = 'genotype', scale.by = 'size', scale.min = 0,
                 cols      = c('lightgrey', 'firebrick'), scale = FALSE,
                 dot.scale = 7) +
    geom_point(aes(size = pct.exp), shape = 21, colour = 'black', stroke = 0.5) +
    coord_flip() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12)) +
    labs(title = 'Top DEGs', x = NULL, y = NULL)
  p_d

  png('plots/Oligodendrocyte_TopDEGs.png',
      width = 4, height = 5, units = 'in', res = 300)
  print(p_d)
  dev.off()
}

# Panel e — rare recDEGs: low pct but significant (expressed in <10% cells)
rare_recDEGs <- WilcoxRNA_sig %>%
  dplyr::filter(pct.1 < 0.1 & pct.2 < 0.1 & pct.1 > 0.03 & pct.2 > 0.03) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n = 10) %>%
  pull(gene)
rare_recDEGs <- intersect(rare_recDEGs, rownames(oligo))

if (length(rare_recDEGs) > 0) {
  p_e <- DotPlot(oligo,
                 features  = rare_recDEGs,
                 group.by  = 'genotype', scale.by = 'size', scale.min = 3,
                 cols      = c('lightgrey', 'firebrick'), scale = FALSE,
                 dot.scale = 7) +
    geom_point(aes(size = pct.exp), shape = 21, colour = 'black', stroke = 0.5) +
    coord_flip() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12)) +
    labs(title = 'Rare DEGs', x = NULL, y = NULL)
  p_e

  png('plots/Oligodendrocyte_RareDEGs.png',
      width = 4, height = 5, units = 'in', res = 300)
  print(p_e)
  dev.off()
}

#==============================================================
# Subclustering (needed for panels g, h, m, n)
#==============================================================
set.seed(42)
oligo <- FindVariableFeatures(oligo, nfeatures = 2000)
oligo <- ScaleData(oligo)
oligo <- RunPCA(oligo, npcs = 30)
ElbowPlot(oligo, ndims = 30)

oligo <- RunUMAP(oligo, dims = 1:10)
oligo <- FindNeighbors(oligo, dims = 1:10)
oligo <- FindClusters(oligo, resolution = 0.1, algorithm = 'leiden',
                      cluster.name = 'oligo_subcluster', random.seed = 1)

#==============================================================
# Save subclustered object
#==============================================================
saveRDS(oligo, 'data/Oligodendrocyte_subclustered.rds')
message('Oligodendrocyte analysis complete.')

#==============================================================
# Panel g — UMAP coloured by subcluster
#==============================================================
oligo <- readRDS('data/Oligodendrocyte_subclustered.rds')

p_g <- DimPlot(oligo, group.by = 'oligo_subcluster', alpha = 0.5) +
  coord_fixed() +
  ggtitle(NULL)
p_g

png('plots/Oligodendrocyte_Subcluster_UMAP.png',
    width = 5, height = 4, units = 'in', res = 300)
print(p_g)
dev.off()

# Also save genotype overlay
p_g2 <- DimPlot(oligo, group.by = 'genotype',
                cols = colmap_genotype, alpha = 0.5) +
  coord_fixed() +
  ggtitle(NULL)
png('plots/Oligodendrocyte_Genotype_UMAP.png',
    width = 5, height = 4, units = 'in', res = 300)
print(p_g2)
dev.off()

#==============================================================
# Panel h — Stacked bar: proportion of each subcluster per genotype
#==============================================================
p_h <- dittoBarPlot(oligo, var = 'oligo_subcluster', group.by = 'genotype', retain.factor.levels = T,
                    color.panel = hue_pal()(4))+
  geom_col(color = 'black') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(title=NULL, x=NULL)
p_h
png('plots/Oligodendrocyte_Composition_Bar.png',
    width = 3, height = 4, units = 'in', res = 300)
print(p_h)
dev.off()

# Save proportion table
wb_prop <- createWorkbook()
addWorksheet(wb_prop, 'SubclusterProportions', gridLines = FALSE)
writeDataTable(wb_prop, 'SubclusterProportions', as.data.frame(table(oligo$sample,oligo$oligo_subcluster,oligo$genotype)))
saveWorkbook(wb_prop, 'results/Oligodendrocyte_SubclusterProportions.xlsx', overwrite = TRUE)

# Bar plots by replicate
countDF <- table(oligo$oligo_subcluster, oligo$genotype, oligo$sample)
countDF <- reshape2::melt(countDF)
colnames(countDF) <- c('oligo_subcluster','genotype','sample','count')
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
              map_signif_level = T, margin_top = 5, vjust = 1.5)+
  facet_wrap(~oligo_subcluster, ncol=5, scales='free')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = colmap_genotype)+
  xlab(NULL)
p1

png('plots/Oligodendrocyte_Composition_Bar_Signif.png', width = 4, height = 3, units = 'in', res = 300)
p1
dev.off()

#==============================================================
# Subcluster marker genes (needed for panels m and n)
#==============================================================
Idents(oligo) <- 'oligo_subcluster'
oligo_markers <- FindAllMarkers(oligo, assay = 'RNA', only.pos = TRUE,
                                logfc.threshold = 1, min.pct = 0.1) %>%
  dplyr::filter(p_val_adj < 0.05)

top10 <- oligo_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(oligo_markers, 'results/Oligodendrocyte_subcluster_markers.csv', row.names = FALSE)

#==============================================================
# Panel m — Heatmap: marker genes for all Subclusters
# Each panel: normalised mean expression (z-score) per genotype × sample
#==============================================================

# Build average expression matrix grouped by genotype × sample × subcluster
avg_expr <- AverageExpression(
  oligo, layer = 'data',
  assays        = 'RNA',
  group.by      = c('genotype','sample','oligo_subcluster'),
  return.seurat = FALSE
)$RNA[unique(top10$gene),]

# Parse column metadata
col_meta <- data.frame(
  col_id   = colnames(avg_expr),
  genotype = substring(sapply(strsplit(colnames(avg_expr), '_'), `[`, 1), 2),
  sample   = sapply(strsplit(colnames(avg_expr), '_'), `[`, 2),
  cluster  = sapply(strsplit(colnames(avg_expr), '_'), `[`, 3),
  stringsAsFactors = FALSE
)

# Sort: cluster → genotype → sample
col_meta <- col_meta[order(as.integer(col_meta$cluster),
                           col_meta$genotype,
                           col_meta$sample), ]
avg_expr <- avg_expr[, col_meta$col_id]

# Row order: genes grouped by source cluster
gene_order <- unique(top10$gene[order(as.integer(as.character(top10$cluster)))])
gene_order <- intersect(gene_order, rownames(avg_expr))
avg_expr   <- avg_expr[gene_order, ]

# z-score rows
avg_scaled <- t(scale(t(avg_expr)))
avg_scaled[avg_scaled >  3] <-  3
avg_scaled[avg_scaled < -3] <- -3

# Colour palettes
n_clusters <- length(unique(col_meta$cluster))
clust_cols <- setNames(hue_pal()(n_clusters), sort(unique(col_meta$cluster)))
geno_cols  <- colmap_genotype

#==============================================================
# Combined marker heatmap across ALL subclusters (full overview)
#==============================================================
row_cluster <- top10$cluster[match(gene_order, top10$gene)]

ha_col_all <- HeatmapAnnotation(
  Genotype = col_meta$genotype,
  col      = list(Genotype = geno_cols),
  annotation_name_side = 'right'
)

ra_row_all <- rowAnnotation(
  Cluster = as.character(row_cluster),
  col     = list(Cluster = clust_cols),
  show_annotation_name = FALSE
)

ht_all <- Heatmap(
  avg_scaled,
  name               = 'z-score',
  col                = circlize::colorRamp2(c(-2, 0, 2),
                                            c('dodgerblue3', 'white', 'firebrick3')),
  top_annotation     = ha_col_all,
  left_annotation    = ra_row_all,
  cluster_rows       = FALSE,
  cluster_columns    = FALSE,
  column_split       = col_meta$cluster,
  row_split          = as.character(row_cluster),
  row_title_gp       = gpar(fontsize = 8, fontface = 'bold'),
  column_title_gp    = gpar(fontsize = 8, fontface = 'bold'),
  row_names_gp       = gpar(fontsize = 8),
  column_names_gp    = gpar(fontsize = 6),
  column_names_rot   = 45,
  show_column_names  = FALSE,
  use_raster         = TRUE
)

png('plots/Oligodendrocyte_Subcluster_MarkerHeatmap_All.png',
    width = 5, height = 2, units = 'in', res = 300)
draw(ht_all)
dev.off()

#==============================================================
# Save subclustered object
#==============================================================
saveRDS(oligo, 'data/Oligodendrocyte_subclustered.rds')
message('Oligodendrocyte analysis complete.')

