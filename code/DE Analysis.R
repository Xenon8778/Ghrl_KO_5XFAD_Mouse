library(dplyr)
library(ggplot2)
library(openxlsx)
library(Seurat)
library(stringr)
library(enrichR)
library(EnhancedVolcano)

#----------------
# Load Dataset
#----------------
so <- readRDS('data/Merged_CellBender_Ghrelin_annotated.rds')
so$ghrelin <- factor(so$ghrelin, levels = c('WT','Ghrl KO'))
so$genotype <- so$ghrelin
so$genotype <- recode_factor(so$genotype,
                                  'WT' = '5XFAD', 'Ghrl KO' = '5XFAD; Ghrl -/-')

table(so$celltype, so$genotype)

ident1 <- "5XFAD; Ghrl -/-"
ident2 <- "5XFAD"

colmap_genotype <- c('5XFAD' = 'steelblue1', '5XFAD; Ghrl -/-' = 'goldenrod1')

dir.create('results/DE_Genotype', showWarnings = FALSE)
dir.create('plots/Volcano_Genotype', showWarnings = FALSE)

#----------------
# Helper: Enrichr (mouse genes, no ortholog conversion needed)
#----------------
setEnrichrSite('Enrichr')

run_enrichr_mouse <- function(glist, wb, sheet_name) {
  dbs <- c('GO_Biological_Process_2025', 'KEGG_2026')

  glist <- as.character(unlist(glist))
  glist <- glist[!is.na(glist) & nchar(glist) > 0]
  if (length(glist) == 0) return(NULL)

  # Capitalize first letter for mouse genes (e.g. Apoe -> APOE for EnrichR)
  glist_upper <- toupper(glist)

  enrichrRes <- enrichr(glist_upper, dbs)

  gobp <- enrichrRes[['GO_Biological_Process_2025']] %>%
    dplyr::filter(P.value < 0.05) %>%
    dplyr::mutate(Term = str_split_i(Term, ' \\(', 1)) %>%
    dplyr::mutate(size = lengths(str_split(Genes, ';'))) %>%
    dplyr::filter(size >= 3)

  kegg <- enrichrRes[['KEGG_2026']] %>%
    dplyr::filter(P.value < 0.05) %>%
    dplyr::mutate(size = lengths(str_split(Genes, ';'))) %>%
    dplyr::filter(size >= 3)

  if (nrow(gobp) > 0) {
    addWorksheet(wb, sheetName = paste0(sheet_name, '_GOBP'), gridLines = FALSE)
    writeDataTable(wb, sheet = paste0(sheet_name, '_GOBP'), x = gobp)
  }
  if (nrow(kegg) > 0) {
    addWorksheet(wb, sheetName = paste0(sheet_name, '_KEGG'), gridLines = FALSE)
    writeDataTable(wb, sheet = paste0(sheet_name, '_KEGG'), x = kegg)
  }

  return(list(gobp = gobp, kegg = kegg))
}

#----------------
# DE Loop: per cell type
#----------------
for (CT in levels(so$celltype)) {
  so1 <- subset(so, subset = celltype == CT)
  n1  <- sum(so1$genotype == ident1)
  n2  <- sum(so1$genotype == ident2)
  message('\n--- ', CT, ' | ', ident1, ': ', n1, ' | ', ident2, ': ', n2, ' ---')

  if (n1 < 10 | n2 < 10) {
    message('Skipping ', CT, ': insufficient cells.')
    next
  }

  wb <- createWorkbook()

  tryCatch({
    # ── Wilcoxon DE ──────────────────────────────────────────────────────────
    WilcoxRNA <- FindMarkers(
      so1,
      group.by       = 'genotype',
      ident.1        = ident1,
      ident.2        = ident2,
      logfc.threshold = 0,       # keep all for volcano; filter below
      min.pct        = 0.05,
      test.use       = 'wilcox'
    ) %>% mutate(gene = rownames(.))

    WilcoxRNA <- WilcoxRNA %>% dplyr::filter(!grepl('^mt-', gene, ignore.case = TRUE)) # remove mitochondrial genes
    WilcoxRNA_sig  <- WilcoxRNA %>% dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 1)
    WilcoxRNA_up   <- WilcoxRNA_sig %>% dplyr::filter(avg_log2FC >  1)
    WilcoxRNA_down <- WilcoxRNA_sig %>% dplyr::filter(avg_log2FC < -1)

    message('  DE sig genes: ', nrow(WilcoxRNA_sig),
            ' (up: ', nrow(WilcoxRNA_up), ', down: ', nrow(WilcoxRNA_down), ')')

    # Write full + split DE tables
    addWorksheet(wb, 'DE_RNA',      gridLines = FALSE)
    writeDataTable(wb, 'DE_RNA',      x = WilcoxRNA)

    addWorksheet(wb, 'DE_RNA_SIG',  gridLines = FALSE)
    writeDataTable(wb, 'DE_RNA_SIG', x = WilcoxRNA_sig)

    addWorksheet(wb, 'DE_RNA_UP',   gridLines = FALSE)
    writeDataTable(wb, 'DE_RNA_UP',  x = WilcoxRNA_up)

    addWorksheet(wb, 'DE_RNA_DOWN', gridLines = FALSE)
    writeDataTable(wb, 'DE_RNA_DOWN', x = WilcoxRNA_down)

    # ── Enrichr ──────────────────────────────────────────────────────────────
    if (nrow(WilcoxRNA_sig) >= 3) {
      run_enrichr_mouse(WilcoxRNA_sig$gene,  wb, 'Enrichr_ALL')
    }
    if (nrow(WilcoxRNA_up) >= 3) {
      run_enrichr_mouse(WilcoxRNA_up$gene,   wb, 'Enrichr_UP')
    }
    if (nrow(WilcoxRNA_down) >= 3) {
      run_enrichr_mouse(WilcoxRNA_down$gene, wb, 'Enrichr_DOWN')
    }

    # ── Save workbook ─────────────────────────────────────────────────────────
    saveWorkbook(wb,
                 paste0('results/DE_Genotype/GhrlKO_vs_WT_', CT, '.xlsx'),
                 overwrite = TRUE)
    message('  Workbook saved.')

    # ── Volcano Plot ──────────────────────────────────────────────────────────
    top_genes <- WilcoxRNA %>%
      dplyr::filter(abs(avg_log2FC) > 1 & p_val_adj < 0.05) %>%
      arrange(p_val_adj) %>%
      pull(gene) %>%
      head(15)

    p_vol <- EnhancedVolcano(
      WilcoxRNA,
      x           = 'avg_log2FC',
      y           = 'p_val_adj',
      lab         = WilcoxRNA$gene,
      selectLab   = top_genes,
      pCutoff     = 0.05,
      FCcutoff    = 1,
      drawConnectors = TRUE,
      col         = c('grey', 'grey', 'steelblue', 'firebrick'),
      title       = paste0(ident1, ' vs ', ident2, ' — ', CT),
      subtitle    = paste0('sig: ', nrow(WilcoxRNA_sig),
                           ' | up: ', nrow(WilcoxRNA_up),
                           ' | down: ', nrow(WilcoxRNA_down))
    )

    png(paste0('plots/Volcano_Genotype/GhrlKO_vs_WT_', CT, '_Volcano.png'),
        width = 8, height = 8, units = 'in', res = 300)
    print(p_vol)
    dev.off()
    message('  Volcano saved.')

  }, error = function(e) {
    cat('ERROR in', CT, ':', conditionMessage(e), '\n')
  })
}

message('\nAll cell types done.')

#----------------
# Summary Heatmap: top 50 DE genes per cell type
#----------------
library(ComplexHeatmap)
library(circlize)
library(scales)

so <- readRDS('data/Merged_CellBender_Ghrelin_annotated.rds')
so$ghrelin <- factor(so$ghrelin, levels = c('WT', 'Ghrl KO'))
so$genotype <- so$ghrelin
so$genotype <- recode_factor(so$genotype,
                             'WT' = '5XFAD', 'Ghrl KO' = '5XFAD; Ghrl -/-')

dir.create('plots/Heatmap_Genotype', showWarnings = FALSE)

# Read DE_RNA_SIG sheet from each saved workbook
xlsx_files <- list.files('results/DE_Genotype', pattern = '^GhrlKO_vs_WT_.+\\.xlsx$', full.names = TRUE)
CTList <- levels(so$celltype)

de_results_list <- sapply(CTList, function(x) {
  fname <- paste0('results/DE_Genotype/GhrlKO_vs_WT_', x, '.xlsx')
  if (read.xlsx(fname, sheet = 'DE_RNA_SIG') %>% nrow() <= 1) {
    message('Skipping ', x, ': DE_RNA_SIG sheet is empty.')
    return(NULL)
  } else{
    read.xlsx(fname, sheet = 'DE_RNA_SIG') %>% mutate(celltype = x)
  }
})
de_results_list <- Filter(Negate(is.null), de_results_list)

for (i in seq_along(de_results_list)) {
  ct <- de_results_list[[i]]$celltype[1]
  DE_df <- de_results_list[[i]]
  message('Cell type: ', ct, ' | DE genes: ', nrow(DE_df))

  # Top 25 up + top 25 down per cell type (= top 50 markers), ordered by p_val_adj
  top_genes_df <- bind_rows(
    bind_rows(DE_df) %>%
      arrange(p_val_adj, .by_group = TRUE) %>%
      dplyr::filter(avg_log2FC >  1) %>% slice_head(n = 25),
    bind_rows(DE_df) %>%
      arrange(p_val_adj, .by_group = TRUE) %>%
      dplyr::filter(avg_log2FC < -1) %>% slice_head(n = 25)
  ) %>% ungroup()

  heatmap_genes <- unique(top_genes_df$gene)
  message('Heatmap genes: ', length(heatmap_genes))

  # Pseudobulk average expression per replicate (celltype x sample)
  so_sub <- subset(so, subset = celltype == ct)
  so_sub <- NormalizeData(so_sub)

  avg_mat <- AggregateExpression(
    so_sub,
    assays   = 'RNA',
    features = heatmap_genes,
    group.by = c('sample','genotype'),
    slot     = 'data'
  )$RNA

  scaled_mat <- ScaleData(avg_mat)  # scale by gene (row) for heatmap

  # Build column metadata by parsing group labels
  col_meta <- data.frame(
    col_id   = colnames(scaled_mat),2,
    sample   = sapply(strsplit(colnames(scaled_mat), '_'), `[`, 1) ,# first field
    genotype = sapply(strsplit(colnames(scaled_mat), '_'), `[`, 2),
    stringsAsFactors = FALSE
  )

  # Order columns: cell type → genotype (5XFAD first) → sample
  genotype_levels <- levels(so$genotype)
  col_meta <- col_meta %>%
    arrange(genotype, sample)
  col_order <- col_meta$col_id
  scaled_mat   <- scaled_mat[, col_order]
  col_meta  <- col_meta[match(col_order, col_meta$col_id), ]

  # Cap extreme values for colour scale
  scaled_mat <- pmin(pmax(scaled_mat, -3), 3)

  # Row annotation: direction and source cell type
  row_annot_df <- top_genes_df %>%
    dplyr::select(gene, avg_log2FC) %>%
    distinct(gene, .keep_all = TRUE) %>%
    mutate(direction = ifelse(avg_log2FC > 0, 'Upregulated', 'Downregulated'))
  rownames(row_annot_df) <- row_annot_df$gene

  # Order rows: by cell type, then up before down
  row_order <- row_annot_df %>%
    mutate(dir_ord = ifelse(direction == 'Upregulated', 1, 2)) %>%
    arrange(dir_ord, gene) %>%
    pull(gene)
  scaled_mat      <- scaled_mat[row_order, ]
  row_annot_df <- row_annot_df[row_order, ]

  # Annotations
  genotype_cols  <- c('5XFAD' = 'steelblue1', '5XFAD; Ghrl -/-' = 'goldenrod1')
  direction_cols <- c('Upregulated' = 'firebrick3', 'Downregulated' = 'steelblue3')

  col_ha <- HeatmapAnnotation(
    Genotype  = col_meta$genotype,
    col = list(Genotype  = genotype_cols),
    annotation_name_side = 'right',
    show_legend = c(TRUE), show_annotation_name = FALSE
  )

  row_ha <- rowAnnotation(
    Direction = row_annot_df$direction,
    col = list(Direction = direction_cols),
    show_legend = FALSE, show_annotation_name = FALSE
  )

  col_fun <- colorRamp2(c(-2.5, 0, 2.5), c('steelblue3', 'white', 'firebrick3'))

  ht <- Heatmap(
    scaled_mat,
    name              = 'Scaled\nexpression',
    col               = col_fun,
    top_annotation    = col_ha,
    left_annotation   = row_ha,
    column_labels = col_meta$sample,
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    column_names_gp   = gpar(fontsize = 10),
    row_names_gp      = gpar(fontsize = 10),
    column_split      = col_meta$genotype,
    row_split         = row_annot_df$direction,
    row_title_gp      = gpar(fontsize = 10, fontface = 'bold'),
    column_title_gp   = gpar(fontsize = 12, fontface = 'bold'),
    border            = TRUE
  )
  draw(ht, merge_legend = TRUE)

  plot_height <- min(10, length(heatmap_genes) * 0.25 + 2)

  png(paste0('plots/Heatmap_Genotype/GhrlKO_vs_WT_Top50_Heatmap_',ct,'.png'),
      width  = 5,
      height = plot_height,
      units  = 'in', res = 300)
  draw(ht, merge_legend = TRUE)
  dev.off()
  message('Summary heatmap saved.')
}


#----------------
# Pathway Dot Plot: Enrichr_ALL_GOBP across cell types
#----------------
library(ComplexHeatmap)
library(circlize)
library(scales)
library(patchwork)

so <- readRDS('data/Merged_CellBender_Ghrelin_annotated.rds')
so$ghrelin <- factor(so$ghrelin, levels = c('WT', 'Ghrl KO'))
so$genotype <- so$ghrelin
so$genotype <- recode_factor(so$genotype,
                             'WT' = '5XFAD', 'Ghrl KO' = '5XFAD; Ghrl -/-')

dir.create('plots/Pathway_Genotype', showWarnings = FALSE)

# Read DE_RNA_SIG sheet from each saved workbook
xlsx_files <- list.files('results/DE_Genotype', pattern = '^GhrlKO_vs_WT_.+\\.xlsx$', full.names = TRUE)

# Read Enrichr_ALL_GOBP sheet from each workbook + join DE direction info
pathway_list <- lapply(xlsx_files, function(f) {
  CT <- sub('.*GhrlKO_vs_WT_(.+)\\.xlsx$', '\\1', basename(f))
  tryCatch({
    gobp <- read.xlsx(f, sheet = 'Enrichr_ALL_GOBP') %>%
      mutate(celltype = CT)

    # Pull avg_log2FC per term: mean FC of the overlapping genes from DE_RNA_SIG
    de   <- read.xlsx(f, sheet = 'DE_RNA_SIG')
    gobp <- gobp %>%
      rowwise() %>%
      mutate(
        overlap_genes = list(toupper(trimws(strsplit(Genes, ';')[[1]]))),
        mean_log2FC   = mean(de$avg_log2FC[toupper(de$gene) %in% overlap_genes],
                             na.rm = TRUE)
      ) %>%
      ungroup() %>%
      dplyr::select(-overlap_genes)

    gobp
  }, error = function(e) { message('Skipping ', CT, ': ', conditionMessage(e)); NULL })
})
pathway_list <- Filter(Negate(is.null), pathway_list)

if (length(pathway_list) > 0) {

  pathway_df <- bind_rows(pathway_list) %>%
    mutate(celltype = factor(celltype, levels = levels(so$celltype)))

  top_n_terms <- 10

  # Build one plot per cell type
  ct_plots  <- list()
  ct_heights <- c()

  for (CT in levels(so$celltype)) {
    df_ct <- pathway_df %>%
      dplyr::filter(celltype == CT) %>%
      arrange(P.value) %>%
      slice_head(n = top_n_terms) %>%
      mutate(Term_label = str_wrap(str_to_title(Term), width = 50),
             Term_label = factor(Term_label, levels = Term_label[order(Odds.Ratio)]))

    if (nrow(df_ct) == 0) next

    p <- ggplot(df_ct, aes(x = Odds.Ratio, y = Term_label,
                           size = size, fill = P.value)) +
      geom_point(shape = 21, color='black') +
      scale_fill_gradient( limits = c(0, 0.05),
        low = 'firebrick', high = 'white',
        name    = 'Adj.\nP-value'
      ) +
      scale_x_log10() +
      # xlim(0,65)+
      scale_size_continuous(limits = c(3, 10), range = c(2,7), name = 'Gene\ncount') +
      theme_classic() +
      theme(
        axis.text.y      = element_text(colour = 'black', size = 10),
        axis.text.x      = element_text(colour = 'black', size = 8),
        legend.position  = 'right'
      ) +
      labs(x = 'Odds Ratio', y = NULL, title = CT)

    ct_plots[[CT]]  <- p
    ct_heights[CT]  <- nrow(df_ct)
  }

  # Stack with patchwork, heights proportional to number of terms, one shared legend
  combined <- wrap_plots(ct_plots, ncol = 1, heights = ct_heights) +
    plot_layout(guides = 'collect') +
    plot_annotation(title = 'GO Biological Process — Ghrl KO vs WT',
                    theme = theme(plot.title = element_text(face = 'bold', size = 12)))
  combined
  total_height <- sum(ct_heights) * 0.3 + length(ct_plots) * 0.8

  png('plots/Pathway_Genotype/GhrlKO_vs_WT_GOBP_Dotplot.png',
      width  = 6,
      height = max(6, total_height),
      units  = 'in', res = 300)
  print(combined)
  dev.off()
  message('Pathway dot plot saved.')

}
