# Aggregate fragment counts per cell type over the pre-IDR peak set in two pseudo-replicates (adapted from ATAC_frags_for_IDR.R)
# changes: pseudo-replicates defined by the 'idr_rep' column of samples.tsv (paper: exp024+exp031+exp042 vs the rest,
#          selected on orig.ident), cell types missing from one replicate skipped, paths from arguments
library(Seurat)
library(Signac)
library(Hmisc)
library(rtracklayer)
source("scripts/pipeline_utils.R")

args = commandArgs(trailingOnly=TRUE)
combined_seurat = readRDS(args[1])
macs_peaks = import.bed(args[2])
out_dir = args[3]
dir.create(out_dir, recursive=TRUE, showWarnings = FALSE)

samples = read_samples()
rep1_ids = samples$sample_id[samples$idr_rep == "1"]
rep2_ids = samples$sample_id[samples$idr_rep == "2"]
if (length(rep1_ids) == 0 || length(rep2_ids) == 0) stop("samples.tsv needs samples with idr_rep 1 and idr_rep 2")

macs_peaks_matrix <- FeatureMatrix(combined_seurat[["WNN"]]@fragments, sep = c("-", "-"), features=macs_peaks, cells = colnames(combined_seurat))
combined_seurat[["macs_peaks"]] <- CreateChromatinAssay(counts = macs_peaks_matrix, sep = c(":", "-"), fragments = combined_seurat[["WNN"]]@fragments)

combined_seurat_1 = subset(combined_seurat, cells = colnames(combined_seurat)[combined_seurat$exp %in% rep1_ids])
combined_seurat_2 = subset(combined_seurat, cells = colnames(combined_seurat)[combined_seurat$exp %in% rep2_ids])

aggregate_rep = function(obj) {
  agg = AggregateExpression(obj, assays = "macs_peaks", group.by = "cell_type", slot = "counts")[["macs_peaks"]]
  as.matrix(agg)
}
combined_seurat_1_agg_atac = aggregate_rep(combined_seurat_1)
combined_seurat_2_agg_atac = aggregate_rep(combined_seurat_2)

shared_types = intersect(colnames(combined_seurat_1_agg_atac), colnames(combined_seurat_2_agg_atac))
missing = setdiff(union(colnames(combined_seurat_1_agg_atac), colnames(combined_seurat_2_agg_atac)), shared_types)
if (length(missing) > 0) warning("cell types present in only one pseudo-replicate (no IDR): ", paste(missing, collapse = ", "))

peak_coords = as.data.frame(matrix(unlist(strsplit(row.names(combined_seurat_1_agg_atac), split = "-")), ncol = 3, byrow = T))
peak_coords$V2 = as.numeric(peak_coords$V2) - 1
peak_coords$name = paste0("embryo_AS_", seq_len(nrow(peak_coords)))

write_rep = function(values, file) {
  rep_df = data.frame(chr = peak_coords$V1, start = peak_coords$V2, end = peak_coords$V3, name = peak_coords$name,
                      score = values, strand = ".", signalValue = values, pValue = values, qValue = values, peak = 100)
  write.table(rep_df, file = file, append = F, quote = F, sep = "\t", row.names = F, col.names = F)
}

for (cell_type in shared_types){
  write_rep(combined_seurat_1_agg_atac[row.names(combined_seurat_1_agg_atac), cell_type], paste0(out_dir, "/", cell_type, "_rep_1.bed"))
  write_rep(combined_seurat_2_agg_atac[row.names(combined_seurat_1_agg_atac), cell_type], paste0(out_dir, "/", cell_type, "_rep_2.bed"))
}

write.table(shared_types, file = paste0(out_dir, "/cell_types.txt"), append = F, quote = F, sep = "\t", row.names = F, col.names = F)

all_types = sort(unique(combined_seurat$cell_type))
k = cbind(rep1 = as.integer(table(factor(combined_seurat_1$cell_type, levels = all_types))),
          rep2 = as.integer(table(factor(combined_seurat_2$cell_type, levels = all_types))))
row.names(k) = all_types
write.table(k, file = paste0(out_dir, "/cell_number_by_type.txt"), append = F, quote = F, sep = "\t", row.names = T, col.names = T)
