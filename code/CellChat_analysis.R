library(Seurat)
library(dplyr)
library(CellChat)
library(patchwork)
library(openxlsx)

colmap_genotype <- c('5XFAD' = 'steelblue1', '5XFAD; Ghrl -/-' = 'goldenrod1')
dir.create('plots/CellChat', showWarnings = FALSE)
dir.create('plots/CellChat/Pathways', showWarnings = FALSE)
dir.create('results/CellChat', showWarnings = FALSE)

#--------------
# Load dataset
#--------------
so <- readRDS('data/Merged_CellBender_Ghrelin_annotated.rds')
so$ghrelin  <- factor(so$ghrelin, levels = c('WT', 'Ghrl KO'))
so$genotype <- recode_factor(so$ghrelin,
                             'WT'      = '5XFAD',
                             'Ghrl KO' = '5XFAD; Ghrl -/-')

ident1 <- '5XFAD'          # WT
ident2 <- '5XFAD; Ghrl -/-' # KO

DefaultAssay(so) <- 'RNA'

CellChatDB <- CellChatDB.mouse

#--------------
# Build CellChat object per genotype
#--------------
build_cellchat <- function(seurat_obj) {
  cc <- createCellChat(object = seurat_obj,
                       group.by = 'celltype',
                       assay = 'RNA')
  cc@DB <- CellChatDB
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = 'triMean')
  cc <- filterCommunication(cc, min.cells = 30)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)
  cc
}

so1 <- subset(so, subset = genotype == ident1)
so2 <- subset(so, subset = genotype == ident2)

cc1 <- build_cellchat(so1)   # WT
cc2 <- build_cellchat(so2)   # KO

saveRDS(cc1, 'results/CellChat/cellchat_WT.rds')
saveRDS(cc2, 'results/CellChat/cellchat_KO.rds')

#--------------
# Merge for comparison
#--------------
object.list <- list(WT = cc1, KO = cc2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat, 'results/CellChat/cellchat_merged.rds')

#--------------
# Compare total interactions/strength
#--------------
p1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
p2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = 'weight')
p1 + p2
png('plots/CellChat/Compare_TotalInteractions.png', width = 6, height = 4, units = 'in', res = 300)
p1 + p2
dev.off()

#--------------
# Differential interaction network
#--------------
png('plots/CellChat/DiffInteraction_Number.png', width = 7, height = 7, units = 'in', res = 300)
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
dev.off()

png('plots/CellChat/DiffInteraction_Weight.png', width = 7, height = 7, units = 'in', res = 300)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = 'weight')
dev.off()

png('plots/CellChat/DiffInteraction_Heatmap_Number.png', width = 6, height = 5, units = 'in', res = 300)
netVisual_heatmap(cellchat)
dev.off()

png('plots/CellChat/DiffInteraction_Heatmap_Weight.png', width = 6, height = 5, units = 'in', res = 300)
netVisual_heatmap(cellchat, measure = 'weight')
dev.off()

#--------------
# Per-group interaction circle plots
#--------------
weight.max <- getMaxWeight(object.list, attribute = c('idents', 'count'))
png('plots/CellChat/CircleNet_perGroup.png', width = 10, height = 6, units = 'in', res = 300)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = TRUE,
                   label.edge = FALSE, edge.weight.max = weight.max[2],
                   edge.width.max = 12, title.name = names(object.list)[i])
}
dev.off()

#--------------
# Information flow comparison (rankNet)
#--------------
p_rank1 <- rankNet(cellchat, mode = 'comparison', stacked = TRUE,  do.stat = TRUE)
p_rank2 <- rankNet(cellchat, mode = 'comparison', stacked = FALSE, do.stat = TRUE)
p_rank1 + p_rank2
png('plots/CellChat/RankNet_InfoFlow.png', width = 8, height = 8, units = 'in', res = 300)
p_rank1 + p_rank2
dev.off()

#--------------
# Signaling role / centrality shifts by cell type
#--------------
num.link <- sapply(object.list, function(x) rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count))
weight.MinMax <- c(min(unlist(num.link)), max(unlist(num.link)))
gg_list <- list()
for (i in 1:length(object.list)) {
  gg_list[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i],
                                                    weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(gg_list, ncol = 2)
png('plots/CellChat/SignalingRole_Scatter.png', width = 10, height = 5, units = 'in', res = 300)
patchwork::wrap_plots(gg_list, ncol = 2)
dev.off()

#--------------
# Shared pathway heatmap comparison
#--------------
pathway.union <- union(object.list[['WT']]@netP$pathways, object.list[['KO']]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(object.list[['WT']], pattern = 'outgoing',
                                         signaling = pathway.union, title = 'WT', width = 7, height = 20)
ht2 <- netAnalysis_signalingRole_heatmap(object.list[['KO']], pattern = 'outgoing',
                                         signaling = pathway.union, title = 'KO', width = 7, height = 20)
png('plots/CellChat/SignalingRole_Heatmap_Outgoing.png', width = 10, height = 10, units = 'in', res = 300)
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, 'cm'))
dev.off()

ht3 <- netAnalysis_signalingRole_heatmap(object.list[['WT']], pattern = 'incoming',
                                         signaling = pathway.union, title = 'WT', width = 7, height = 20,
                                         color.heatmap = 'GnBu')
ht4 <- netAnalysis_signalingRole_heatmap(object.list[['KO']], pattern = 'incoming',
                                         signaling = pathway.union, title = 'KO', width = 7, height = 20,
                                         color.heatmap = 'GnBu')
png('plots/CellChat/SignalingRole_Heatmap_Incoming.png', width = 10, height = 10, units = 'in', res = 300)
ComplexHeatmap::draw(ht3 + ht4, ht_gap = unit(0.5, 'cm'))
dev.off()

#--------------
# Per-pathway differential circle plots (union of pathways)
#--------------
for (pw in pathway.union) {
  tryCatch({
    pw.max <- getMaxWeight(object.list, slot.name = 'netP', attribute = pw)
    png(paste0('plots/CellChat/Pathways/Pathway_', gsub('[^A-Za-z0-9]', '_', pw), '.png'),
        width = 10, height = 6, units = 'in', res = 300)
    par(mfrow = c(1, 2), xpd = TRUE)
    for (i in 1:length(object.list)) {
      netVisual_aggregate(object.list[[i]], signaling = pw, layout = 'circle',
                          edge.weight.max = pw.max, signaling.name = paste(pw, names(object.list)[i]))
    }
    dev.off()
  }, error = function(e) cat('ERROR pathway', pw, ':', conditionMessage(e), '\n'))
}

#--------------
# Identify significantly differential pathways (rankNet stats)
#--------------
rank_df <- rankNet(cellchat, mode = 'comparison', stacked = TRUE, do.stat = TRUE,
                   return.data = TRUE)$signaling.contribution
sig_pathways <- rank_df %>% dplyr::filter(pvalues < 0.05) %>%
  dplyr::pull(name) %>% unique()
write.csv(rank_df, 'results/CellChat/Pathway_DifferentialStats.csv', row.names = FALSE)

#--------------
# Per-pathway differential circle plots (significant pathways only)
#--------------
dir.create('plots/CellChat/DiffPathways', showWarnings = FALSE)
for (pw in sig_pathways) {
  tryCatch({
    pw.max <- getMaxWeight(object.list, slot.name = 'netP', attribute = pw)
    png(paste0('plots/CellChat/DiffPathways/DiffPathway_', gsub('[^A-Za-z0-9]', '_', pw), '.png'),
        width = 10, height = 6, units = 'in', res = 300)
    par(mfrow = c(1, 2), xpd = TRUE)
    for (i in 1:length(object.list)) {
      netVisual_aggregate(object.list[[i]], signaling = pw, layout = 'circle',
                          edge.weight.max = pw.max, signaling.name = paste(pw, names(object.list)[i]))
    }
    dev.off()
  }, error = function(e) cat('ERROR pathway', pw, ':', conditionMessage(e), '\n'))
}

#--------------
# Ligand-receptor level differential expression (identifyOverExpressed by group)
#--------------
cellchat_all <- createCellChat(object = so, group.by = 'celltype', assay = 'RNA')
cellchat_all@DB <- CellChatDB
cellchat_all <- subsetData(cellchat_all)
cellchat_all <- addMeta(cellchat_all, meta = so@meta.data)
cellchat_all <- setIdent(cellchat_all, ident.use = 'celltype')
cellchat_all@meta$datasets <- so$genotype

wb <- createWorkbook()
for (CT in levels(so$celltype)) {
  tryCatch({
    res <- identifyOverExpressedGenes(cellchat_all,
                                      idents.use = CT, group.dataset = 'datasets',
                                      pos.dataset = ident2, features.name = paste0(CT, '_up'),
                                      only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1,
                                      thresh.p = 0.05, do.DE = TRUE)
    net.diff <- res@var.features[[paste0(CT, '_up', '.info')]]
    if (!is.null(net.diff) && nrow(net.diff) > 0) {
      addWorksheet(wb, sheetName = substr(CT, 1, 31), gridLines = FALSE)
      writeDataTable(wb, sheet = substr(CT, 1, 31), x = net.diff)
    }
  }, error = function(e) cat('ERROR DE', CT, ':', conditionMessage(e), '\n'))
}
saveWorkbook(wb, 'results/CellChat/LR_DE_byCelltype.xlsx', overwrite = TRUE)

message('CellChat differential analysis complete.')
