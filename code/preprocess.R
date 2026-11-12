library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(SeuratWrappers)
library(dittoSeq)
library(scales)
library(ggsignif)

colmap_ghrelin <- c('WT' = 'steelblue1', 'Ghrl KO' = 'goldenrod1')
colmap_genotype <- c('5XFAD' = 'steelblue1', '5XFAD; Ghrl -/-' = 'goldenrod1')

#--------------
# Custom CellBender H5 reader
#--------------
library(hdf5r)
f <- H5File$new('data/cellbender_outs/1-FG109_cb_outs_filtered.h5', mode = 'r')
f$ls(recursive = TRUE)
f$close()

library(hdf5r)
library(Matrix)
read_cellbender_h5 <- function(path) {
  f <- H5File$new(path, mode = 'r')

  counts   <- f[['matrix/data']][]
  indices  <- f[['matrix/indices']][]
  indptr   <- f[['matrix/indptr']][]
  barcodes <- f[['matrix/barcodes']][]
  features <- f[['matrix/features/name']][]
  dims     <- f[['matrix/shape']][]

  f$close()

  mat <- sparseMatrix(
    i = indices + 1,
    p = indptr,
    x = as.numeric(counts),
    dims = dims,
    dimnames = list(make.unique(features), barcodes)  # fix duplicate gene names
  )
  return(mat)
}

#--------------
# Load data
#--------------

# 1-FG109 | WT | 14mo
so_FG109 <- CreateSeuratObject(counts = read_cellbender_h5('data/cellbender_outs/1-FG109_cb_outs_filtered.h5'), assay = "RNA")
so_FG109$sample   <- "FG-109"
so_FG109$ghrelin  <- "WT"
so_FG109$age   <- 14
so_FG109

# 3-FG156 | WT | 14mo
so_FG156 <- CreateSeuratObject(counts = read_cellbender_h5('data/cellbender_outs/3-FG156_cb_outs_filtered.h5'), assay = "RNA")
so_FG156$sample   <- "FG-156"
so_FG156$ghrelin  <- "WT"
so_FG156$age   <- 14
so_FG156

# 5-FG176 | WT | 13mo
so_FG176 <- CreateSeuratObject(counts = read_cellbender_h5('data/cellbender_outs/5-FG176_cb_outs_filtered.h5'), assay = "RNA")
so_FG176$sample   <- "FG-176"
so_FG176$ghrelin  <- "WT"
so_FG176$age   <- 13
so_FG176

# 2-FG143 | Ghrl KO | 14mo
so_FG143 <- CreateSeuratObject(counts = read_cellbender_h5('data/cellbender_outs/2-FG143_cb_outs_filtered.h5'), assay = "RNA")
so_FG143$sample   <- "FG-143"
so_FG143$ghrelin  <- "Ghrl KO"
so_FG143$age   <- 14
so_FG143

# 4-FG174 | Ghrl KO | 13mo
so_FG174 <- CreateSeuratObject(counts = read_cellbender_h5('data/cellbender_outs/4-FG174_cb_outs_filtered.h5'), assay = "RNA")
so_FG174$sample   <- "FG-173"
so_FG174$ghrelin  <- "Ghrl KO"
so_FG174$age   <- 13
so_FG174

# 6-FG175 | Ghrl KO | 13mo
so_FG175 <- CreateSeuratObject(counts = read_cellbender_h5('data/cellbender_outs/6-FG175_cb_outs_filtered.h5'), assay = "RNA")
so_FG175$sample   <- "FG-175"
so_FG175$ghrelin  <- "Ghrl KO"
so_FG175$age   <- 13
so_FG175

#--------------
# Merge
#--------------
combinedSO <- merge(
  x = so_FG109,
  y = list(so_FG156, so_FG176, so_FG143, so_FG174, so_FG175),
  add.cell.ids = c("FG109", "FG156", "FG176", "FG143", "FG174", "FG175")
)
combinedSO

saveRDS(combinedSO, 'data/Merged_CellBender_Ghrelin.rds')

#--------------
# Filtering
#--------------
so <- readRDS('data/Merged_CellBender_Ghrelin.rds')
table(so$ghrelin, so$sample)
so

# Quality Control
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")

p1 <- VlnPlot(
  object = so,
  features = c("nCount_RNA","nFeature_RNA","percent.mt"),
  ncol = 3,
  pt.size = 0
)
p1
png('plots/QC_violins.png', width = 6, height = 4, units = 'in', res = 300)
p1
dev.off()

soFilt <- subset(
  x = so,
  subset = percent.mt < 5 &
    nCount_RNA < 50000 &
    nCount_RNA > 1000 &
    nFeature_RNA > 500
)
soFilt

library(future)
plan("multisession", workers = 4)
options(future.globals.maxSize= 1.6e+10)
plan("sequential")

# Normalization, PCA and Harmony
# soFilt <- JoinLayers(soFilt)
soFilt <- NormalizeData(soFilt)
soFilt <- FindVariableFeatures(soFilt)
soFilt <- ScaleData(soFilt)
soFilt <- RunPCA(soFilt)
ElbowPlot(soFilt, reduction = 'pca', ndims = 50)

soFilt <- IntegrateLayers(
  object = soFilt, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE
)
soFilt <- JoinLayers(soFilt)

# UMAP
soFilt <- FindNeighbors(soFilt, reduction = 'harmony', dims = 1:30,
                        return.neighbor = T)
soFilt <- RunUMAP(soFilt, umap.method = 'uwot2',
                  nn.name = 'RNA.nn')

DimPlot(soFilt, reduction = 'umap', group.by = 'sample')
DimPlot(soFilt, reduction = 'umap', group.by = 'ghrelin')

saveRDS(soFilt, 'data/Merged_CellBender_Ghrelin_UMAP.rds')

#--------------
# Annotation
#--------------
soUMAP <- readRDS('data/Merged_CellBender_Ghrelin_UMAP.rds')
DimPlot(soUMAP, group.by = 'ghrelin', reduction = "umap")

# Clustering
set.seed(1)
soUMAP <- FindNeighbors(soUMAP, reduction = 'harmony', dims = 1:30)
soUMAP <- FindClusters(soUMAP, algorithm = 'leiden', resolution = 1.2,
                       cluster.name = 'leiden', random.seed = 1)

DimPlot(soUMAP, group.by = 'leiden',  label = TRUE, reduction = "umap") + NoLegend()

# Find Marker genes
markers <- FindAllMarkers(soUMAP, group.by = 'leiden', assay = 'RNA',
                          only.pos = T, logfc.threshold = 1, min.pct = 0.1)
write.csv(markers %>% filter(p_val_adj < 0.01),'results/leiden_cluster_markers.csv')

markers <- read.csv('results/leiden_cluster_markers.csv')
top20 <- markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, -avg_log2FC) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup()

clipr::write_clip(markers %>% dplyr::filter(cluster == 10) %>%
                    head(100) %>% pull(gene))

# Read
DimPlot(soUMAP, group.by = 'ghrelin',  label = TRUE, repel = TRUE, reduction = "umap")
DimPlot(soUMAP, group.by = 'leiden',  label = TRUE, repel = TRUE, reduction = "umap",
        raster = T) + NoLegend()

feats <- c('Nes','Sox2','Pax6', # Stem
           'Stmn2','Rbfox3','Map2', # Neuron
           'Foxg1', 'Slc17a7','Slc17a6', # Excitatory
           'Gad1', 'Gad2', 'Pvalb',# GABA
           'Chat', 'Slc5a7', # Cholinergic
           'Slc1a3', 'S100b','Aqp4','Aldh1l1','Gfap',  #Astro 2
           'Itgam','Cd68','P2ry12', 'Tmem119', # Microglia
           'Mbp','Plp1', 'Mog', 'Mag','Cnp', # Oligo
           'Pdgfra', # OPC,
           'Foxj1', 'Pifo','Dynlrb2','Trp73','Tmem212', # Ependymal
           'Mpz','Sox10', # Schwann
           'Ptprc', # Immune
           'Clic6','Ttr','Kl', # Choroid plexus (CP)
           'Pecam1', # Endothelial
           'Col3a1','Dcn','Fn1' # Fibroblast
)

DotPlot(soUMAP, group.by = 'leiden', features = feats, dot.scale = 5,
        cols = 'RdBu', cluster.idents = T, scale.max = 80)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  coord_flip()

soUMAP$celltype = soUMAP$leiden
soUMAP$celltype = recode_factor(soUMAP$celltype,
                                '1' = 'Oligodendrocyte',
                                '2' = 'GABAergic',
                                '3' = 'GABAergic',
                                '4' = 'Oligodendrocyte',
                                '5' = 'GABAergic',
                                '6' = 'Glutaminergic',
                                '7' = 'Microglia',
                                '8' = 'Glutaminergic',
                                '9' = 'Glutaminergic',
                                '10' = 'Glutaminergic',
                                '11' = 'Glutaminergic',
                                '12' = 'Glutaminergic',
                                '13' = 'GABAergic',
                                '14' = 'Astrocyte',
                                '15' = 'OPC',
                                '16' = 'Microglia',
                                '17' = 'Glutaminergic',
                                '18' = 'Glutaminergic',
                                '19' = 'GABAergic',
                                '20' = 'GABAergic',
                                '21' = 'Microglia',
                                '22' = 'GABAergic',
                                '23' = 'Glutaminergic',
                                '24' = 'GABAergic',
                                '25' = 'GABAergic',
                                '26' = 'Fibroblast',
                                '27' = 'CP',
                                '28' = 'GABAergic',
                                '29' = 'Glutaminergic',
                                '30' = 'GABAergic',
                                '31' = 'GABAergic',
                                '32' = 'Astrocyte',
                                '33' = 'Glutaminergic',
                                '34' = 'Glutaminergic',
                                '35' = 'Glutaminergic',
                                '36' = 'OPC',
                                '37' = 'Glutaminergic',
                                '38' = 'GABAergic',
                                '39' = 'Ependymal',
                                '40' = 'GABAergic',
                                '41' = 'Oligodendrocyte',
                                '42' = 'Cholinergic',
                                '43' = 'GABAergic',
                                '44' = '.')
DimPlot(soUMAP, group.by = 'celltype') + coord_fixed()+
  ggtitle('Cell Types')

soUMAP$celltype <- factor(soUMAP$celltype,
                          levels = c('GABAergic','Glutaminergic','Cholinergic',
                                     'Astrocyte','Microglia','Oligodendrocyte','OPC',
                                     'Ependymal','CP','Fibroblast'))
DimPlot(soUMAP, group.by = 'celltype') + coord_fixed()+
  ggtitle('Cell Types')


table(soUMAP$celltype, soUMAP$ghrelin)

feats <- c('Rbfox3','Map2', # Neuron
           'Gad1', 'Gad2', # GABA
           'Foxg1', 'Slc17a7','Slc17a6', # Glutaminergic
           'Chat', 'Slc5a7', # Cholinergic
           'Slc1a3', 'S100b','Aqp4','Aldh1l1','Gfap',  #Astro
           'Itgam','Cd68','P2ry12', 'Tmem119', # Microglia
           'Mbp','Plp1', 'Mog', 'Mag','Cnp', # Oligo
           'Pdgfra', # OPC,
           'Foxj1', 'Pifo','Dynlrb2','Tmem212', # Ependymal
           'Clic6','Ttr','Kl', # Choroid plexus (CP)
           'Col3a1','Dcn','Fn1' # Fibroblast
)
p1 <- DotPlot(soUMAP, group.by = 'celltype', features = feats, dot.scale = 5,
              cols = 'RdBu',scale.max = 50)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(x=NULL,y=NULL)
p1

saveRDS(soUMAP, 'data/Merged_CellBender_Ghrelin_annotated.rds')


#--------------
# Plots
#--------------
soAnnot <- readRDS('data/Merged_CellBender_Ghrelin_annotated.rds')
soAnnot$ghrelin <- factor(soAnnot$ghrelin, levels = c('WT','Ghrl KO'))
soAnnot$genotype <- soAnnot$ghrelin
soAnnot$genotype <- recode_factor(soAnnot$genotype,
                                  'WT' = '5XFAD', 'Ghrl KO' = '5XFAD; Ghrl -/-')

countDF <- as.data.frame.matrix(table(soAnnot$celltype, soAnnot$ghrelin))
countDF
percDF <- t(countDF) / rowSums(t(countDF)) * 100
percDF


wb <- createWorkbook()
addWorksheet(wb, sheetName = 'CellNumbers', gridLines = FALSE)
writeDataTable(wb = wb, sheet = 'CellNumbers', x = countDF, rowNames = T)
addWorksheet(wb, sheetName = 'CellPercentages', gridLines = FALSE)
writeDataTable(wb = wb, sheet = 'CellPercentages', x = percDF, rowNames = T)
saveWorkbook(wb, paste0('results/CellNumbersPercentages.xlsx'),
             overwrite = TRUE)

# Leiden UMAP
p1 <- DimPlot(soAnnot, group.by = 'leiden', alpha = 0.2, pt.size = 0.05) +
  coord_fixed() +
  ggtitle('Leiden Clusters')
p1
png('plots/Leiden_UMAP.png', width = 7, height = 6, units = 'in', res = 300)
p1
dev.off()

# Markers used for annotation
feats <- c('Nes','Sox2','Pax6', # Stem
           'Stmn2','Rbfox3','Map2', # Neuron
           'Foxg1', 'Slc17a7','Slc17a6', # Excitatory
           'Gad1', 'Gad2', 'Pvalb',# GABA
           'Chat', 'Slc5a7', # Cholinergic
           'Slc1a3', 'S100b','Aqp4','Aldh1l1','Gfap',  #Astro 2
           'Itgam','Cd68','P2ry12', 'Tmem119', # Microglia
           'Mbp','Plp1', 'Mog', 'Mag','Cnp', # Oligo
           'Pdgfra', # OPC,
           'Foxj1', 'Pifo','Dynlrb2','Trp73','Tmem212', # Ependymal
           'Mpz','Sox10', # Schwann
           'Ptprc', # Immune
           'Clic6','Ttr','Kl', # Choroid plexus (CP)
           'Pecam1', # Endothelial
           'Col3a1','Dcn','Fn1' # Fibroblast
)

p1 <- DotPlot(soAnnot, group.by = 'leiden', features = feats, dot.scale = 5,
              cols = 'RdBu', cluster.idents = T)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(x=NULL,y=NULL)
p1

png('plots/Cluster_marker_dotplot.png', width = 10, height = 10, units = 'in', res = 300)
p1
dev.off()


# Celltype UMAP
p1 <- DimPlot(soAnnot, group.by = 'celltype', alpha = 0.2, pt.size = 0.05) +
  coord_fixed() +
  ggtitle('Cell Types')
p1
png('plots/Celltype_UMAP.png', width = 7, height = 6, units = 'in', res = 300)
p1
dev.off()

# Genotype UMAP
p1 <- DimPlot(soAnnot, group.by = 'genotype', alpha = 0.2, pt.size = 0.05,
              cols = colmap_genotype) +
  coord_fixed() +
  ggtitle('Genotypes')
p1
png('plots/Genotype_UMAP.png', width = 7, height = 6, units = 'in', res = 300)
p1
dev.off()

# Celltype UMAP by Genotype
p1 <- DimPlot(soAnnot, group.by = 'celltype', alpha = 0.2, pt.size = 0.05,
              split.by = 'genotype') +
  coord_fixed() +
  ggtitle('Cell Types by Genotype')
p1
png('plots/Celltype_by_Genotype_UMAP.png', width = 11, height = 6, units = 'in', res = 300)
p1
dev.off()

# Celltype marker genes
FindMarkers(soAnnot, ident.1 = 'Glutaminergic', group.by = 'celltype',
            only.pos=T)

feats <- c('Rbfox3', # Neuron
           'Gad1', 'Gad2', # GABA
           'Slc17a7','Satb2','Camk2a', # Glutaminergic
           'Chat', 'Slc5a7', # Cholinergic
           'Slc1a3', 'S100b','Aqp4','Aldh1l1','Gfap',  #Astro
           'Itgam','Cd68','P2ry12', 'Tmem119', # Microglia
           'Mbp','Plp1', 'Mog', 'Mag','Cnp', # Oligo
           'Pdgfra', # OPC,
           'Foxj1', 'Pifo','Dynlrb2','Tmem212', # Ependymal
           'Clic6','Ttr','Kl', # Choroid plexus (CP)
           'Col3a1','Dcn','Fn1' # Fibroblast
)
p1 <- DotPlot(soAnnot, group.by = 'celltype', features = feats, dot.scale = 5,
              cols = 'RdBu',scale.max = 50)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(x=NULL,y=NULL)
p1

png('plots/Celltype_marker_dotplot.png', width = 5, height = 9, units = 'in', res = 300)
p1
dev.off()

# Bar plots
table(soAnnot$celltype, soAnnot$genotype)
p1 <- dittoBarPlot(soAnnot, var = 'celltype', group.by = 'genotype', retain.factor.levels = T,
                   color.panel = hue_pal()(10))+
  geom_col(color = 'black') + NoLegend() + ggtitle(NULL)+
  xlab('Genotype')
p2 <- dittoBarPlot(soAnnot, var = 'celltype', group.by = 'genotype', retain.factor.levels = T,
                   color.panel = hue_pal()(10), scale = 'count')+
  geom_col(color = 'black') + ggtitle(NULL)+
  xlab('Genotype')
p1+p2
png('plots/Genotype_celltype_barplot.png', width = 5, height = 6, units = 'in', res = 300)
p1+p2
dev.off()

# Bar plots by replicate
countDF <- table(soAnnot$celltype, soAnnot$genotype, soAnnot$sample)
countDF <- reshape2::melt(countDF)
colnames(countDF) <- c('celltype','genotype','sample','count')
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
              map_signif_level = T,margin_top = 0.8, vjust = 1.5)+
  facet_wrap(~celltype, ncol=5, scales='free')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = colmap_genotype)+
  xlab(NULL)
p1

png('plots/Genotype_celltype_sample_barplot_signif.png', width = 8, height = 6, units = 'in', res = 300)
p1
dev.off()

percDFsorted <- percDF
percDFsorted$celltype <- factor(percDFsorted$celltype, levels = unique(sort_by(percDFsorted$celltype,-percDFsorted$Percentage)))
p1 <- percDFsorted %>% ggplot(aes(x=celltype, y=Percentage, fill=genotype, group = genotype))+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2,
               position = position_dodge(0.8))+
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.6,
               position = position_dodge(0.8))+
  geom_jitter(aes(shape = genotype),
              position = position_dodge(0.8))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = colmap_genotype)+
  xlab(NULL)
p1
png('plots/Genotype_celltype_sample_barplot.png', width = 8, height = 3, units = 'in', res = 300)
p1
dev.off()
