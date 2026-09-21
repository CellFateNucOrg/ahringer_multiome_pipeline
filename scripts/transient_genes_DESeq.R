# Transient (pre-gastrulation) genes: pseudo-bulk DESeq2 of late-cell-cycle nuclei of pre- vs post-gastrulation lineages
# (adapted from transient_genes_DESeq.20241113.R)
# changes: cell-cycle calls read from the object, lineages absent from the data are dropped, heatmap column order
#          no longer relies on the paper's cluster numbers, WormCAT / commented code removed
library(Seurat)
library(Signac)
library(Hmisc)
library(RColorBrewer)
library(ComplexHeatmap)
library(DESeq2)
library(ggplot2)
source("scripts/pipeline_utils.R")

plot_dir_w_date = make_plot_dir("transient_genes")
dir.create("transient_genes", showWarnings = FALSE)

args = commandArgs(trailingOnly = TRUE)
combined_seurat = readRDS(args[1])
wb_id_common_names = read.table(paste0(canon_prefix(), ".filtered.gene_name.txt"))
combined_seurat$cell_cycle = combined_seurat$cell_cycle_rf

# lineages outside gastrulation, labelled pre-gastrulation (early) or post-gastrulation (late)
lineage_stage = c("ABx" = "early", "ABxx" = "early", "ABxxx" = "early", "ABxxxxx" = "late", "ABxxxxxx" = "late",
                  "EMS" = "early", "E" = "early", "Exx" = "late", "Exxx" = "late",
                  "MS" = "early", "MSx" = "early", "MSxxx" = "late", "MSxxxx" = "late",
                  "C" = "early", "Cx" = "early", "Cxxx" = "late",
                  "Dx" = "late", "DxxCxxx" = "late")

# sample up to 200 late cell-cycle nuclei per lineage
sample_cell_n = 200
set.seed(100)
cell_lineage_barcodes_all = list()
for (cell_l in names(lineage_stage)) {
  bc = colnames(combined_seurat)[combined_seurat$cell_lineage == cell_l & combined_seurat$cell_cycle == "late"]
  if (length(bc) == 0) next
  cell_lineage_barcodes_all[[cell_l]] = if (length(bc) <= sample_cell_n) bc else sample(bc, size = sample_cell_n, replace = F)
}
lineages = names(cell_lineage_barcodes_all)
if (sum(lineage_stage[lineages] == "early") < 2 || sum(lineage_stage[lineages] == "late") < 2) {
  stop("need at least 2 pre- and 2 post-gastrulation lineages with late cell-cycle nuclei")
}
message("lineages used: ", paste(lineages, collapse = ", "))
combined_seurat_sample = subset(combined_seurat, cells = unlist(cell_lineage_barcodes_all))

# pseudo-bulk UMI counts per lineage
countData_lineage = as.matrix(AggregateExpression(combined_seurat_sample, assays = "RNA", group.by = "cell_lineage", slot = "counts")[["RNA"]])
countData_lineage = round(countData_lineage[, lineages])

colData <- data.frame(row.names = colnames(countData_lineage),
                      cell_lineage = lineages,
                      pregastr_late = relevel(factor(unname(lineage_stage[lineages])), ref = "late"))

set.seed(42)
dds <- DESeqDataSetFromMatrix(countData = countData_lineage, colData = colData, design = ~pregastr_late)
dds <- DESeq(dds)
res <- results(dds, alpha = 0.001, altHypothesis = "greater", lfcThreshold = 1)
write.table(as.data.frame(res), file = "transient_genes/DESeq2_pregastrulation_vs_postgastrulation.txt", quote = F, sep = "\t")

# transient genes with high / low expression in post-gastrulation lineages
k = counts(dds, normalized = TRUE)
late_cols = colData$pregastr_late == "late"
high_late = row.names(k[which(apply(k[, late_cols, drop = FALSE], 1, function(x){quantile(x, 0.25)}) > 25),])

DE_id_lfc_2_se7 = row.names(res[res$padj < 0.001 & !is.na(res$padj) & res$log2FoldChange > 2 & res$lfcSE < 0.7,])
DE_id_lfc_2_se7_hl = DE_id_lfc_2_se7[DE_id_lfc_2_se7 %in% high_late]
DE_id_lfc_2_se7_ll = DE_id_lfc_2_se7[DE_id_lfc_2_se7 %nin% high_late]
message(length(DE_id_lfc_2_se7), " transient genes")

write.table(wb_id_common_names[wb_id_common_names$V2 %in% DE_id_lfc_2_se7, 1], file = "transient_genes/transient_genes.all.txt", quote = F, row.names = F, col.names = "gene_name")
write.table(wb_id_common_names[wb_id_common_names$V2 %in% DE_id_lfc_2_se7_hl, 1], file = "transient_genes/transient_genes.high_postgastrulation_exp.txt", quote = F, row.names = F, col.names = "gene_name")
write.table(wb_id_common_names[wb_id_common_names$V2 %in% DE_id_lfc_2_se7_ll, 1], file = "transient_genes/transient_genes.low_postgastrulation_exp.txt", quote = F, row.names = F, col.names = "gene_name")

# average expression heatmaps
cell_lineage_levels <- c("twocell", "ABx", "ABxx", "ABxxx", "ABxxxx", "ABxxxxx", "ABxxxxxx",
                         "EMS", "EandMS", "E", "Ex", "Exx", "Exxx", "MS", "MSx", "MSxx", "MSxxx", "MSxxxx",
                         "C", "Cx", "Cxx", "Cxxx", "D", "Dx", "DxxCxxx", "unassigned", "P2", "P3", "P4", "Z2Z3")
stage_levels <- c("twocell", "4cell", "8cell", "15cell", "26cell", "46cell", "87morecell", "Pcell", "Z2Z3", "unassigned")

avg_heatmap = function(group, levels = NULL) {
  m = as.matrix(AverageExpression(combined_seurat, assays = "RNA", group.by = group, features = DE_id_lfc_2_se7)[["RNA"]])
  ord = if (is.null(levels)) sort(colnames(m)) else c(intersect(levels, colnames(m)), setdiff(sort(colnames(m)), levels))
  m = m[, ord, drop = FALSE]
  draw(Heatmap(log10(m + 1), name = "transient gene exp", cluster_columns = F, cluster_rows = T, use_raster = F,
               row_dend_reorder = TRUE, column_title_rot = 90, show_row_dend = F, show_row_names = F,
               column_title = group), padding = unit(c(8, 5, 10, 4), "mm"))
}

if (length(DE_id_lfc_2_se7) > 1) {
  pdf(file = paste0(plot_dir_w_date, "/transient_genes_heatmap.pdf"), width = 16, height = 16)
  avg_heatmap("cell_lineage", cell_lineage_levels)
  avg_heatmap("cell_type")
  avg_heatmap("stage", stage_levels)
  dev.off()
}
