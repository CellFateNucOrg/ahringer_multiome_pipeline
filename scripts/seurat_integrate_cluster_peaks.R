# Add embryo stage metadata and the final IDR peak set (+ FRiP) to the annotated object
# (adapted from seurat_integrate_cluster_peaks.20240913.R; unused colour vectors removed, fragment paths read from the object)
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(glmGamPoi)
library(Signac)
library(rtracklayer)
library(Hmisc)
source("scripts/pipeline_utils.R")

set.seed(42)
args = commandArgs(trailingOnly=TRUE)

combined_seurat_rds = args[1]
merged_peaks_bed = args[2]
new_ATAC_assay = args[3]
final_seurat_rds = args[4]
gene_annotation = args[5]

combined_seurat <- readRDS(combined_seurat_rds)
ct = unique(combined_seurat$cell_type)

# embryo stage per cell type
cell_types_by_stage =list("twocell" = "twocell",
                          "4cell" = c("ABx", "earlyEMS", "EMS"),
                          "8cell" = c("ABax", "ABpx", "earlyEandMS", "E", "MS", "earlyC", "C"),
                          "15cell" = c("earlyABxxx", "ABxxx", "ABala", "ABalp", "ABara", "ABarp", "ABpxa", "ABpxp", "earlyEx", "earlyMSx", "MSx", "MSa", "MSp", "earlyCx"),
                          "26cell" = c("earlyABxxxx", "ABxxxx", "ABalaa", "ABalap", "ABalpa", "ABalpp", "ABaraa", "ABarap", "ABarpa", "ABarpp", "ABpxaa", "ABpxap", "ABpxpa", "ABpxpp", "Ex", "earlyMSxx", "Ca", "Cp", "D"),
                          "46cell" = c("MSxx", "MSxa", "MSxp", "earlyExx", "earlyCxx", "Cxa", "Cxp",
                                       ct[intersect(grep(pattern = "ABxxxxx", ct), grep(pattern = "ABxxxxxx", ct, invert = T))]),
                          "87morecell" = c(ct[grep(pattern = "ABxxxxxx", ct)],
                                           "Exx", "Exxx",
                                           ct[grep(pattern = "MSxxx", ct)],
                                           ct[grep(pattern = "MSxxxx", ct)],
                                           ct[intersect(grep(pattern = "Cxxx", ct), grep(pattern = "DxxCxxx", ct, invert = T))],
                                           ct[grep(pattern = "DxxCxxx", ct)],
                                           "Dx"),
                          "Pcell" = c("P2", "P3", "P4"),
                          "Z2Z3" = "Z2Z3",
                          "unassigned" = ct[grep(pattern = "unassigned", ct)]
)

stage_metadata = combined_seurat$cell_type
for (stage in names(cell_types_by_stage)) {
  stage_metadata[combined_seurat$cell_type %in% cell_types_by_stage[[stage]]] = stage
}
combined_seurat$stage = stage_metadata

# add final peak set
merged_peaks = import.bed(merged_peaks_bed)
merged_peaks_matrix_filtered <- FeatureMatrix(combined_seurat[["WNN"]]@fragments, sep = c("-", "-"), features=merged_peaks, cells = colnames(combined_seurat))
combined_seurat[[new_ATAC_assay]] <- CreateChromatinAssay(
  counts = merged_peaks_matrix_filtered,
  sep = c(":", "-"),
  fragments = combined_seurat[["WNN"]]@fragments)

DefaultAssay(combined_seurat) <- new_ATAC_assay
combined_seurat <- RunTFIDF(combined_seurat)

# FRiP
combined_seurat$tot_frag = total_fragments(combined_seurat, new_ATAC_assay)
combined_seurat = FRiP(object = combined_seurat, assay = new_ATAC_assay, total.fragments = "tot_frag", col.name = "FRiP")

saveRDS(combined_seurat, file = final_seurat_rds)
write.table(combined_seurat[[]], file = sub("\\.rds$", ".table.txt", final_seurat_rds), quote = F, sep = "\t", row.names = T, col.names = T)
