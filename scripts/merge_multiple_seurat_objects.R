# Merge per-sample objects, re-cluster on RNA (SCT), add ATAC counts over the bulk peak set from sample-prefixed
# fragment files, and write per-cluster fragments for MACS2 (adapted from merge_multiple_seurat_objects.20240430.R)
# changes: any number of samples (comma-separated list) instead of exactly 7; unused ATAC h5 loading removed
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(glmGamPoi)
library(Signac)
library(rtracklayer)
source("scripts/pipeline_utils.R")

set.seed(42)
args = commandArgs(trailingOnly=TRUE)

exp_ids = strsplit(args[1], ",")[[1]]
merged_macs = args[2]
fragments_out_dir = args[3]
annotation = args[4]

out_prefix = paste0("seurat_objects/combined_all_samples_SCT.", annotation)
plot_dir = make_plot_dir("merge_samples")
cR_output = "data/sc/cellranger_arc/"

objs = lapply(exp_ids, function(id) {
  o <- readRDS(paste0("seurat_objects/", id, ".", annotation, ".postSoupX.ATAC_MACS.rds"))
  o@meta.data$exp <- id
  o
})

gcdata <- merge_seurat_list(objs, add.cell.ids = exp_ids, merge.data = T)
rm(objs); gc()

VariableFeatures(gcdata[["SCT"]]) <- rownames(gcdata[["SCT"]]@scale.data)
DefaultAssay(gcdata) = "SCT"

gcdata <- RunPCA(gcdata, seed.use = 42, features = VariableFeatures(object = gcdata), npcs = 100)
gcdata <- FindNeighbors(gcdata, dims = 1:50)
gcdata <- FindClusters(gcdata, resolution = 3, random.seed = 0)
gcdata <- RunUMAP(gcdata, seed.use = 42, dims = 1:50, return.model=TRUE)

pdf(file=paste0(plot_dir, "/all_samples_SCT.", annotation, ".umap_split.50PCA.pdf"), width=18, height = 8)
print(ElbowPlot(gcdata, ndims = 100))
print(DimPlot(gcdata, reduction = "umap", group.by = "orig.ident"))
print(DimPlot(gcdata, reduction = "umap", split.by = "orig.ident"))
print(DimPlot(gcdata, reduction = "umap", label = TRUE ))
dev.off()

plot_key_markers(gcdata, "umap", paste0(plot_dir, "/all_samples_SCT.", annotation, ".key_lineage_markers.50PCA.pdf"))

# add ATAC data: counts over the bulk peak set, using fragment files with sample-prefixed barcodes
merged_peaks = import.bed(merged_macs)
fragment_list = lapply(exp_ids, function(id) {
  CreateFragmentObject(paste0(cR_output, id, "/outs/atac_fragments_barcode_corrected.tsv.gz"),
                       cells = colnames(gcdata)[startsWith(colnames(gcdata), paste0(id, "_"))])
})

merged_peaks_matrix_filtered <- FeatureMatrix(fragment_list, sep = c("-", "-"), features=merged_peaks, cells = colnames(gcdata))
merged_peaks_matrix_filtered1 <- merged_peaks_matrix_filtered[,colnames(merged_peaks_matrix_filtered) %in% colnames(gcdata)]

gcdata[["ATAC"]] <- CreateChromatinAssay(
  counts = merged_peaks_matrix_filtered1,
  sep = c(":", "-"),
  fragments = fragment_list)

pdf(file=paste0(plot_dir, "/all_samples_SCT.", annotation, ".all_cells.MACS_peaks.nfeature_atac_rna.pdf"), width=24, height = 8)
p1 = DimPlot(gcdata, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4)
p2 = FeaturePlot(gcdata, features = "nFeature_ATAC", min.cutoff = 0, max.cutoff = 2000, order=T) & scale_colour_gradientn(colours = feature_colours)
p3 = FeaturePlot(gcdata, features = "nFeature_RNA", min.cutoff = 0, max.cutoff = 1000, order=T) & scale_colour_gradientn(colours = feature_colours)
print(p1 + p2 + p3)
p2 = FeaturePlot(gcdata, features = "nCount_ATAC", min.cutoff = 0, max.cutoff = 3000, order=T) & scale_colour_gradientn(colours = feature_colours)
p3 = FeaturePlot(gcdata, features = "nCount_RNA", min.cutoff = 0, max.cutoff = 2000, order=T) & scale_colour_gradientn(colours = feature_colours)
print(p1 + p2 + p3)
dev.off()

saveRDS(gcdata, file = paste0(out_prefix, ".all_cells.MACS_peaks.rds"))
write.table(gcdata[[]], file = paste0(out_prefix, ".all_cells.MACS_peaks.table.txt"), quote = F, sep = "\t", row.names = T, col.names = T)

# fragment files per cluster (for MACS2 / bigWigs)
DefaultAssay(gcdata) <- "ATAC"
dir.create(fragments_out_dir, recursive = T, showWarnings = F)
unlink(file.path(fragments_out_dir, "*.bed"))
SplitFragments(gcdata, group.by = "seurat_clusters", outdir = fragments_out_dir, append = T)

write_cells_per_group(gcdata, "seurat_clusters", paste0(out_prefix, ".all_cells.MACS_peaks.cells_per_cluster.txt"))
