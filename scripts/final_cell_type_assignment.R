# Final cell type / lineage assignment (adapted from final_cell_type_assignment.20240911.R)
# changes:
#  * unannotated barcodes detected as "not in the AB / E_MS_C_D tables" instead of "barcode contains 'exp'"
#  * P-cell sub-clusters ranked by m1m2 accessibility for any number of sub-clusters (paper assumed exactly 3)
#  * late-cluster table optional; cell types relabelled 'unassigned' set by UNASSIGN_CELL_TYPES (paper: Cxx)
#  * output files without date stamps
library(ggplot2)
library(dplyr)
library(Seurat)
library(Signac)
library(Hmisc)
library(stringr)
library(multimode)
library(rtracklayer)

source("scripts/second_derivative_binarisation.R")
source("scripts/pipeline_utils.R")

plot_dir_w_date = make_plot_dir("final_cell_types")

args = commandArgs(trailingOnly=TRUE)
seurat_object_file = args[1]
annotation = args[2]
late_clusters = if (length(args) >= 3) args[3] else ""
unassign_types = strsplit(env_chr("UNASSIGN_CELL_TYPES", ""), "[ ,]+")[[1]]
out_prefix = paste0("seurat_objects/wt_all.all_samples.all_cells.", annotation)

AB_barcodes_file = paste0("cell_type_annotation_all/AB.cell_assignment.", annotation, ".txt")
E_MS_C_D_barcodes_file = paste0("cell_type_annotation_all/E_MS_C_D.cell_assignment.", annotation, ".txt")
P_granule_file = "external_data/P_granule_transcripts.txt"

combined_seurat = readRDS(seurat_object_file)

# cell type information for the somatic early lineages
read_assignment = function(f) if (file.exists(f)) read.table(f) else data.frame(cell_type = character(0))
AB_barcodes = read_assignment(AB_barcodes_file)
E_MS_C_D_barcodes = read_assignment(E_MS_C_D_barcodes_file)

# remove underscores from cell type names, rename "early_E_MS" to "earlyEandMS" (avoid mix-up with "earlyEMS")
AB_barcodes$cell_type <- gsub('_','',AB_barcodes$cell_type)
E_MS_C_D_barcodes$cell_type[E_MS_C_D_barcodes$cell_type == "early_E_MS"] = "earlyEandMS"
E_MS_C_D_barcodes$cell_type <- gsub('_','',E_MS_C_D_barcodes$cell_type)

all_barcodes = as.character(combined_seurat$seurat_clusters)
names(all_barcodes) = colnames(combined_seurat)
ab_idx = match(row.names(AB_barcodes), colnames(combined_seurat))
all_barcodes[ab_idx[!is.na(ab_idx)]] = AB_barcodes$cell_type[!is.na(ab_idx)]
emscd_idx = match(row.names(E_MS_C_D_barcodes), colnames(combined_seurat))
all_barcodes[emscd_idx[!is.na(emscd_idx)]] = E_MS_C_D_barcodes$cell_type[!is.na(emscd_idx)]
# barcodes whose sub-cluster got no assignment (NA) fall back to their cluster number
na_idx = is.na(all_barcodes)
all_barcodes[na_idx] = as.character(combined_seurat$seurat_clusters)[na_idx]
combined_seurat$cell_type = unname(all_barcodes)

# subcluster the putative 2-cell cluster and retain only the cpg-3-enriched subclusters as "twocell"
combined_seurat_avgExp_cpg3 = as.numeric(as.data.frame(AverageExpression(combined_seurat, slot = "data", assays = "RNA", features = "cpg-3")))
putative_twocell = which(combined_seurat_avgExp_cpg3 == max(combined_seurat_avgExp_cpg3))-1

combined_seurat = FindSubCluster(combined_seurat, cluster = putative_twocell, graph.name = "wsnn", resolution = 0.5, subcluster.name = "twocell")
combined_seurat$twocell = as.character(combined_seurat$twocell)
twocell_sucl_names = unique(combined_seurat@meta.data$twocell)[grep(unique(combined_seurat@meta.data$twocell), pattern = "_")]
twocell_subclusters = c()
for (i in twocell_sucl_names){
  if ("cpg-3" %in% row.names(FindMarkers(combined_seurat, ident.1 = i, group.by = "twocell", only.pos = T, assay = "RNA", logfc.threshold = 2))){
    twocell_subclusters = append(twocell_subclusters, i)
  }
}
combined_seurat$cell_type[combined_seurat$twocell %in% twocell_subclusters] = "twocell"

# putative Z2/Z3 cluster
combined_seurat_avgExp_prg2 = as.numeric(as.data.frame(AverageExpression(combined_seurat, slot = "data", assays = "RNA", features = "prg-2")))
putative_z2z3 = which(combined_seurat_avgExp_prg2 == max(combined_seurat_avgExp_prg2))-1
combined_seurat$cell_type[combined_seurat$twocell %in% putative_z2z3] = "Z2Z3"

### P cells: enrichment of the P-granule module, then sub-clustering
p_granule_genes = read.table(P_granule_file)
combined_seurat = AddModuleScore(combined_seurat, features = list(p_granule_genes$V1), assay = "RNA")
combined_seurat$p_granules = combined_seurat$Cluster1
combined_seurat$Cluster1 = NULL

all_clusters = unique(combined_seurat@meta.data$twocell)
all_clusters = all_clusters[all_clusters %nin% twocell_subclusters]
combined_seurat_avg_p_granule_exp = sapply(all_clusters, function(cl) mean(combined_seurat$p_granules[combined_seurat$twocell == cl]))

pdf(file=paste0(plot_dir_w_date, "/P_granule_binarisation.pdf"), width=8, height = 8)
loc <- locmodes(combined_seurat_avg_p_granule_exp, mod0 = 2, display = T)
low_high_exp_loc = second_derivative_binarization(combined_seurat_avg_p_granule_exp, loc$cbw$bw, loc$locations, "P-granules")
dev.off()
p_cell_clusters <- all_clusters[which(combined_seurat_avg_p_granule_exp >= low_high_exp_loc[2])]
message("P-cell clusters: ", paste(p_cell_clusters, collapse = ","))
combined_seurat$cell_type[combined_seurat$twocell %in% p_cell_clusters] = "Pcell"

Idents(combined_seurat) = "cell_type"
combined_seurat = FindSubCluster(combined_seurat, cluster = "Pcell", graph.name = "wsnn", resolution = 0.3, subcluster.name = "pcell")
combined_seurat$cell_type = as.character(combined_seurat$pcell)
pcell_sucl_names = sort(unique(combined_seurat$cell_type[grep("^Pcell_", combined_seurat$cell_type)]))

# accessibility of germline-specific m1m2 promoters distinguishes P2 < P3 < P4
m1m2_elements = import.bed("external_data/reg_elements_all.elegans.gl_specific.m1m2.bed")
m1m2_elements_matrix <- FeatureMatrix(combined_seurat[["WNN"]]@fragments, sep = c("-", "-"), features=m1m2_elements, cells = colnames(combined_seurat))
combined_seurat[["m1m2_elements"]] <- CreateChromatinAssay(counts = m1m2_elements_matrix, sep = c(":", "-"), fragments = combined_seurat[["WNN"]]@fragments)

mean_counts = sapply(pcell_sucl_names, function(cl) mean(combined_seurat$nCount_m1m2_elements[combined_seurat$cell_type == cl]))
ranked = names(sort(mean_counts))
if (length(ranked) != 3) warning(length(ranked), " P-cell sub-clusters found (paper: 3); lowest m1m2 -> P2, highest -> P4, others -> P3")
if (length(ranked) >= 1) combined_seurat$cell_type[combined_seurat$cell_type == ranked[1]] = "P2"
if (length(ranked) >= 2) combined_seurat$cell_type[combined_seurat$cell_type == ranked[length(ranked)]] = "P4"
if (length(ranked) >= 3) combined_seurat$cell_type[combined_seurat$cell_type %in% ranked[2:(length(ranked)-1)]] = "P3"

pdf(file=paste0(plot_dir_w_date, "/Pcell_annotation.pdf"), width=12, height = 6)
par(mfrow=c(1,2))
germ = intersect(c("P2", "P3", "P4", "Z2Z3"), unique(combined_seurat$cell_type))
boxplot(lapply(germ, function(g) combined_seurat$nFeature_m1m2_elements[combined_seurat$cell_type == g]),
        outline = F, main = "Accessible m1m2 sites in germline precursors", las = 1, names = germ)
plot(1, type = "n", xlim = c(0, max(combined_seurat$nFeature_m1m2_elements)), ylim = c(0, 0.1), las = 1,
     main = "Accessible m1m2 sites in germline precursors", xlab = "nFeature m1m2", ylab = "density")
for (k in seq_along(germ)) {
  v = combined_seurat$nFeature_m1m2_elements[combined_seurat$cell_type == germ[k]]
  if (length(v) > 1) lines(density(v, bw = 3), col = k, lwd = 2, lty = k)
}
legend("topright", legend = germ, fill = seq_along(germ), bty = "n")
dev.off()

# late clusters (manual table: cluster number -> lineage label)
cell_type_seurat_cluster = combined_seurat$cell_type
reassigned_cell_type = cell_type_seurat_cluster
if (late_clusters != "" && file.exists(late_clusters)) {
  lines_ok = grep("^\\s*[^#[:space:]]", readLines(late_clusters), value = TRUE)
  if (length(lines_ok) > 0) {
    reassigned_clusters_table = read.table(text = lines_ok, colClasses = "character")
    for (cluster_n in reassigned_clusters_table$V1){
      reassigned_cell_type[cell_type_seurat_cluster == cluster_n] = paste0(reassigned_clusters_table[reassigned_clusters_table$V1 == cluster_n, 2], cluster_n)
    }
  }
}
combined_seurat$cell_type = reassigned_cell_type
cl_n = max(as.numeric(combined_seurat$seurat_clusters))-1
combined_seurat$cell_type[combined_seurat$cell_type %in% as.character(c(0:cl_n))] = "unassigned"

# cell types judged unreliable after inspection (paper: Cxx - high ATAC counts, low RNA complexity)
combined_seurat$cell_type[combined_seurat$cell_type %in% unassign_types] = "unassigned"

# lineage information
cell_lineage = combined_seurat$cell_type
cell_lineage = unlist(strsplit(cell_lineage, "early"))
cell_lineage = cell_lineage[cell_lineage != ""]

new_cell_lineage = cell_lineage
new_cell_lineage[str_detect(cell_lineage, "AB[a-z]{2}$")] = "ABxx"
new_cell_lineage[str_detect(cell_lineage, "AB[a-z]{3}$")] = "ABxxx"
new_cell_lineage[str_detect(cell_lineage, "AB[a-z]{4}$")] = "ABxxxx"
new_cell_lineage[str_detect(cell_lineage, "AB[a-z]{5}")] = "ABxxxxx"
new_cell_lineage[str_detect(cell_lineage, "AB[a-z]{6}")] = "ABxxxxxx"
new_cell_lineage[str_detect(cell_lineage, "MS[a-z]{1}$")] = "MSx"
new_cell_lineage[str_detect(cell_lineage, "MS[a-z]{2}$")] = "MSxx"
new_cell_lineage[str_detect(cell_lineage, "MS[a-z]{3}")] = "MSxxx"
new_cell_lineage[str_detect(cell_lineage, "MS[a-z]{4}")] = "MSxxxx"
new_cell_lineage[str_detect(cell_lineage, "C[a-z]{1}$")] = "Cx"
new_cell_lineage[str_detect(cell_lineage, "C[a-z]{2}$")] = "Cxx"
new_cell_lineage[str_detect(cell_lineage, "Cxxx")] = "Cxxx"
new_cell_lineage[str_detect(cell_lineage, "D[a-z]{1}$")] = "Dx"
new_cell_lineage[str_detect(cell_lineage, "DxxCxxx")] = "DxxCxxx"
new_cell_lineage[str_detect(cell_lineage, "unassigned")] = "unassigned"
combined_seurat$cell_lineage = as.character(new_cell_lineage)

pdf(file=paste0(plot_dir_w_date, "/wt.all_samples.all_cells.", annotation, ".UMAP.pdf"), width=16, height = 8)
p1 = DimPlot(combined_seurat, reduction = "umapWNN", group.by = "cell_type", order = T, label = T, label.size = 4) & NoLegend()
p2 = DimPlot(combined_seurat, reduction = "umapWNN", group.by = "cell_lineage", order = T, label = T, label.size = 4) & NoLegend()
print(p1 | p2)
dev.off()

drop_cols = intersect(c("p_granules", "twocell", "SCT_snn_res.3", "WNN_snn_res.3"), colnames(combined_seurat@meta.data))
combined_seurat@meta.data[drop_cols] <- NULL

# fragment files per cell lineage and cell type (for MACS2 and tracks)
DefaultAssay(combined_seurat) <- "WNN"
for (cell_annotation in c("cell_lineage", "cell_type")) {
  d = paste0("cluster_specific_peaks/round_2_all_samples/", annotation, "/", cell_annotation, "/bed")
  dir.create(d, recursive = T, showWarnings = F)
  unlink(file.path(d, "*.bed"))
  SplitFragments(combined_seurat, group.by = cell_annotation, outdir = d, append = T)
  write_cells_per_group(combined_seurat, cell_annotation, paste0(out_prefix, ".cells_per_", cell_annotation, ".txt"))
}

saveRDS(combined_seurat, file = paste0(out_prefix, ".rds"))
write.table(combined_seurat[[]], file = paste0(out_prefix, ".table.txt"), quote = F, sep = "\t", row.names = T, col.names = T)
print(table(combined_seurat$cell_type))
