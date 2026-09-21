# Add the cluster-specific peak set and cluster on RNA + ATAC with WNN
# (adapted from merge_multiple_seurat_objects.cluster_peaks_step_2.20240430.R; sample ids no longer needed)
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(glmGamPoi)
library(Signac)
library(rtracklayer)
source("scripts/pipeline_utils.R")

set.seed(42)
args = commandArgs(trailingOnly=TRUE)

combined_seurat_rds = args[1]
merged_peaks_bed = args[2]
new_ATAC_assay = args[3]
annotation = args[4]
out_rds = sub("\\.all_cells\\.MACS_peaks\\.rds$", paste0(".", new_ATAC_assay, ".rds"), combined_seurat_rds)
plot_dir_w_date = make_plot_dir("wnn")

combined_seurat <- readRDS(combined_seurat_rds)
merged_peaks = import.bed(merged_peaks_bed)

merged_peaks_matrix_filtered <- FeatureMatrix(combined_seurat[["ATAC"]]@fragments, sep = c("-", "-"), features=merged_peaks, cells = colnames(combined_seurat))

combined_seurat[[new_ATAC_assay]] <- CreateChromatinAssay(
  counts = merged_peaks_matrix_filtered,
  sep = c(":", "-"),
  fragments = combined_seurat[["ATAC"]]@fragments)

# cluster barcodes based on ATAC data
DefaultAssay(combined_seurat) <- new_ATAC_assay
combined_seurat <- RunTFIDF(combined_seurat)
combined_seurat <- FindTopFeatures(combined_seurat, min.cutoff = 'q0')
combined_seurat <- RunSVD(combined_seurat)
combined_seurat <- RunUMAP(combined_seurat, reduction = 'lsi', dims = 2:50, reduction.name = "umapATAC")
combined_seurat <- FindNeighbors(combined_seurat, reduction = 'lsi', dims = 2:50)
combined_seurat <- FindClusters(combined_seurat, verbose = T, algorithm = 3, resolution = 3)

# cluster barcodes based on both modalities
combined_seurat <- FindMultiModalNeighbors(combined_seurat, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
combined_seurat <- RunUMAP(combined_seurat, nn.name = "weighted.nn", reduction.name = "umapWNN")
combined_seurat <- FindClusters(combined_seurat, graph.name = "wsnn", algorithm = 3, verbose = T, resolution = 4)

p1 <- DimPlot(combined_seurat, reduction = "umap", group.by = "wsnn_res.4", label = TRUE, label.size = 4)
p2 <- DimPlot(combined_seurat, reduction = "umapATAC", group.by = "wsnn_res.4", label = TRUE, label.size = 4)
p3 <- DimPlot(combined_seurat, reduction = "umapWNN", group.by = "wsnn_res.4", label = TRUE, label.size = 4)

pdf(file=paste0(plot_dir_w_date, "/all_samples_SCT.", annotation, ".umap_WNN.pdf"), width=24, height = 8)
print(p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
dev.off()

pdf(file=paste0(plot_dir_w_date, "/all_samples_SCT.", annotation, ".umap_split.WNN.pdf"), width=24, height = 8)
p1 = DimPlot(combined_seurat, reduction = "umapWNN", label = TRUE, label.size = 4)
p2 = DimPlot(combined_seurat, reduction = "umapWNN", group.by = "exp")
print(p1 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)) | p2  & theme(plot.title = element_text(hjust = 0.5)))
print(DimPlot(combined_seurat, reduction = "umapWNN", split.by = "exp") & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
dev.off()

DefaultAssay(combined_seurat) = "SCT"
plot_key_markers(combined_seurat, "umapWNN", paste0(plot_dir_w_date, "/all_samples_SCT.", annotation, ".key_lineage_markers.WNN.pdf"))

saveRDS(combined_seurat, file = out_rds)
