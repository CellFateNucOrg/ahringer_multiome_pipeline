# Marker-based annotation of AB and E/MS/C/D sub-clusters (adapted from lineage_specific_annotation.20240207.R)
# changes: plots in plots/lineage_annotation; a marker whose expression cannot be binarised by the
#          second-derivative method gets a simple max-relative fallback (with a warning) instead of stopping
library(Seurat)
library(Signac)
library(ggplot2)
library(multimode)
library(Hmisc)
library(stringr)

source("scripts/second_derivative_binarisation.R")
source("scripts/pipeline_utils.R")

plot_dir_w_date = make_plot_dir("lineage_annotation")

args = commandArgs(trailingOnly=TRUE)
seurat_object_file = args[1]
output_dir = args[2]
annotation = args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

### scoring function
scoring_markers = function(x, y){
  score_cluster_cell = length(x)
  malus_value = score_cluster_cell * 0.1
  kill_value = malus_value * 5
  for (marker_gene in c(1:length(x))){
    if (x[marker_gene] == y[marker_gene]){
      score_cluster_cell = score_cluster_cell
    }
    else if (abs(x[marker_gene] - y[marker_gene]) == 2){
      score_cluster_cell = score_cluster_cell - kill_value
    }
    else if (abs(x[marker_gene] - y[marker_gene]) == 1){
      score_cluster_cell = score_cluster_cell - malus_value
    }
  }
  return(score_cluster_cell)
}

combined_seurat = readRDS(seurat_object_file)

# reduce resolution to facilitate identification of AB, E, MS, C and D lineages based on marker genes
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5, graph.name = "wsnn")

AB_mkrs = c("sea-2","T25E12.6","sre-45")
E_MS_C_D_mkrs = c("F28H6.8","elt-7","tbx-35", "ceh-51", "hnd-1", "nlp-11", "cbs-1")

if (file.exists(paste0(output_dir, "/all_clusters_markers.", annotation, ".txt"))){
  all_mrkrs = read.table(paste0(output_dir, "/all_clusters_markers.", annotation, ".txt"))
} else {
  all_mrkrs = FindAllMarkers(combined_seurat, assay = "SCT", logfc.threshold = 0.5, only.pos = T)
  write.table(all_mrkrs, file = paste0(output_dir, "/all_clusters_markers.", annotation, ".txt"), quote = F, sep = "\t", row.names = T, col.names = T)
}

AB_clusters = unique(all_mrkrs[all_mrkrs$gene %in% AB_mkrs,6])
E_MS_C_D_clusters = unique(all_mrkrs[all_mrkrs$gene %in% E_MS_C_D_mkrs,6])
message("AB clusters: ", paste(AB_clusters, collapse = ","), " | E/MS/C/D clusters: ", paste(E_MS_C_D_clusters, collapse = ","))

pdf(file=paste0(plot_dir_w_date, "/early_lineage_markers.", annotation, ".pdf"), width=18, height = 9)
print(FeaturePlot(combined_seurat, reduction = "umapWNN", order = T, features = AB_mkrs, ncol = 3)  & scale_colour_gradientn(colours = feature_colours))
print(FeaturePlot(combined_seurat, reduction = "umapWNN", order = T, features = E_MS_C_D_mkrs, ncol = 4)  & scale_colour_gradientn(colours = feature_colours))
AB_UMAP = DimPlot(combined_seurat, reduction = "umapWNN", label = T, sizes.highlight = 0.4, pt.size = 0.4, cells.highlight =
          colnames(combined_seurat[,combined_seurat$seurat_clusters %in% AB_clusters])) & NoLegend()
E_MS_C_D_UMAP = DimPlot(combined_seurat, reduction = "umapWNN", label = T, sizes.highlight = 0.4, pt.size = 0.4, cells.highlight =
          colnames(combined_seurat[,combined_seurat$seurat_clusters %in% E_MS_C_D_clusters])) & NoLegend()
print(AB_UMAP | E_MS_C_D_UMAP)
dev.off()

lineage_clusters = list("AB" = AB_clusters, "E_MS_C_D" = E_MS_C_D_clusters)

binarise_gene = function(values, gene_name) {
  tryCatch({
    loc <- locmodes(values, mod0 = 2, display = F)
    second_derivative_binarization(values, loc$cbw$bw, loc$locations, gene_name)
  }, error = function(e) {
    warning("binarisation failed for ", gene_name, " (", conditionMessage(e), "); using 25%/60% of max as thresholds")
    m = max(values)
    if (m <= 0) c(Inf, Inf) else c(0.25 * m, 0.6 * m)
  })
}

for (lineage in c("AB", "E_MS_C_D")){
  if (length(lineage_clusters[[lineage]]) == 0) { warning("no cluster found for lineage ", lineage); next }
  markers_table = read.table(paste0("data/external_data/markers_UMAP_annotation.", lineage, ".20240207.tsv"))
  combined_seurat_subcluster = subset(combined_seurat, idents = lineage_clusters[[lineage]])
  DefaultAssay(combined_seurat_subcluster) = "WNN"
  combined_seurat_subcluster <- FindMultiModalNeighbors(combined_seurat_subcluster, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50), weighted.nn.name = "clusterwnn")
  combined_seurat_subcluster <- FindClusters(combined_seurat_subcluster, resolution = 5.5, graph.name = "wsnn")
  combined_seurat_subcluster <- RunUMAP(combined_seurat_subcluster, nn.name = "clusterwnn", reduction.name = "clusterumapWNN")
  pdf(file = paste0(plot_dir_w_date, "/", lineage, "_subclustering.", annotation, ".pdf"), width = 22, height = 6, useDingbats = F)
  p1 = DimPlot(combined_seurat_subcluster, reduction = "clusterumapWNN", label = T) & NoLegend()
  p2 = DimPlot(combined_seurat_subcluster, reduction = "clusterumapWNN", label = T, group.by = "wsnn_res.4") & NoLegend()
  p3 = DimPlot(combined_seurat, reduction = "umapWNN", label = T, group.by = "wsnn_res.4") & NoLegend()
  print(p1 | p2 | p3)
  dev.off()
  # identify modes of expression and divide genes into high/low/no expression
  cl_n = max(as.numeric(combined_seurat_subcluster$seurat_clusters))-1
  markers_table = markers_table[row.names(markers_table) %in% row.names(combined_seurat_subcluster[["RNA"]]), , drop = FALSE]
  combined_seurat_subcluster_avgExp = as.data.frame(AverageExpression(combined_seurat_subcluster, slot = "data", assays = "RNA", features = row.names(markers_table)))
  combined_seurat_subcluster_avgExp = combined_seurat_subcluster_avgExp[row.names(markers_table), , drop = FALSE]
  combined_seurat_subcluster_avgExp_log = log10(combined_seurat_subcluster_avgExp + 1)
  combined_seurat_subcluster_avgExp_quant = combined_seurat_subcluster_avgExp_log
  pdf(file = paste0(plot_dir_w_date, "/", lineage, "_marker_binarisation.", annotation, ".pdf"), width = 8, height = 8, useDingbats = F)
  for (gene in c(1:dim(combined_seurat_subcluster_avgExp_log)[1])){
    values = as.numeric(combined_seurat_subcluster_avgExp_log[gene,])
    low_high_exp_loc = binarise_gene(values, row.names(combined_seurat_subcluster_avgExp_log)[gene])
    combined_seurat_subcluster_avgExp_quant[gene,which(values < low_high_exp_loc[1])] <- 0
    combined_seurat_subcluster_avgExp_quant[gene,which(values >= low_high_exp_loc[1] & values < low_high_exp_loc[2])] <- 1
    combined_seurat_subcluster_avgExp_quant[gene,which(values >= low_high_exp_loc[2])] <- 2
  }
  dev.off()
  # binarised expression as metadata
  for (gene in c(1:dim(combined_seurat_subcluster_avgExp_quant)[1])){
    marker_gene_exp = as.numeric(combined_seurat_subcluster$seurat_clusters)-1
    for (cluster_n in c(0:cl_n)){
      marker_gene_exp[as.numeric(combined_seurat_subcluster$seurat_clusters)-1 == cluster_n] = combined_seurat_subcluster_avgExp_quant[gene, cluster_n+1]
    }
    combined_seurat_subcluster@meta.data[[row.names(combined_seurat_subcluster_avgExp_quant)[gene]]] = marker_gene_exp
  }
  pdf(file = paste0(plot_dir_w_date, "/", lineage, "_marker_binarisation.UMAP.", annotation, ".pdf"), width = 22, height = 6, useDingbats = F)
  # default assay stays WNN so that 'gene' resolves to the binarised metadata column and 'rna_gene' to expression
  for (gene in row.names(markers_table)){
    p1 = FeaturePlot(combined_seurat_subcluster, reduction = "clusterumapWNN", order = T, features = paste0("rna_", gene)) & scale_colour_gradientn(colours = c("grey90", "cyan", "deepskyblue", "darkorange", "darkviolet"))
    p2 = FeaturePlot(combined_seurat_subcluster, reduction = "clusterumapWNN", order = T, features = gene) & scale_colour_gradientn(colours = c("grey90", "cyan", "deepskyblue", "darkorange", "darkviolet"))
    p3 = DimPlot(combined_seurat_subcluster, reduction = "clusterumapWNN", label = T) & NoLegend()
    print(p1 | p2 | p3)
  }
  dev.off()
  # score subclusters based on markers expression
  scoring_cluster = matrix(NA, nrow = dim(markers_table)[2], ncol = dim(combined_seurat_subcluster_avgExp_quant)[2], dimnames = list(colnames(markers_table), colnames(combined_seurat_subcluster_avgExp_quant)))
  for (cluster_name in c(0:cl_n)){
    for (cell_type in c(1:dim(markers_table)[2])){
      scoring_cluster[cell_type, cluster_name+1] = scoring_markers(markers_table[,cell_type],combined_seurat_subcluster_avgExp_quant[, cluster_name+1])
    }
  }
  write.table(scoring_cluster, file = paste0(output_dir, "/", lineage, "_clusters_scoring.", annotation, ".txt"), quote = F, row.names = T, col.names = T)
  # assign subclusters to cell types
  n_markers = dim(markers_table)[1]
  marker_assignment = rep(NA, cl_n + 1)
  for (cluster_n in c(1:dim(scoring_cluster)[2])){
    candidate_cell = which(scoring_cluster[,cluster_n] == max(scoring_cluster[,cluster_n]))
    if (max(scoring_cluster[,cluster_n]) >= n_markers*0.7 & length(candidate_cell) == 1){
      marker_assignment[cluster_n] = names(candidate_cell)
    }
    else if (scoring_cluster[,cluster_n][order(scoring_cluster[,cluster_n])][dim(scoring_cluster)[1]] > 0 &
             scoring_cluster[,cluster_n][order(scoring_cluster[,cluster_n])][dim(scoring_cluster)[1]] - scoring_cluster[,cluster_n][order(scoring_cluster[,cluster_n])][dim(scoring_cluster)[1]-1] > n_markers*0.5){
      marker_assignment[cluster_n] = names(candidate_cell)
    }
  }
  new_metadata = as.numeric(combined_seurat_subcluster$seurat_clusters)
  for (cluster_n in c(1:length(marker_assignment))){
    if (!is.na(marker_assignment[cluster_n])){
      new_metadata[as.numeric(combined_seurat_subcluster$seurat_clusters) == cluster_n] = marker_assignment[cluster_n]
    }
    else{
      new_metadata[as.numeric(combined_seurat_subcluster$seurat_clusters) == cluster_n] = NA
    }
  }
  combined_seurat_subcluster$cell_assignment = new_metadata
  p1 = DimPlot(combined_seurat_subcluster, reduction = "clusterumapWNN", order = T, label = T, group.by = "cell_assignment")
  p2 = DimPlot(combined_seurat_subcluster, reduction = "clusterumapWNN", order = T, label = T, group.by = "wsnn_res.5.5") & NoLegend()
  p3 = DimPlot(combined_seurat_subcluster, reduction = "clusterumapWNN", order = T, label = T, group.by = "wsnn_res.4") & NoLegend()
  pdf(file = paste0(plot_dir_w_date, "/", lineage, "_full_annotation.before_rescue.UMAP.", annotation, ".pdf"), width = 14, height = 6, useDingbats = F)
  print(p1 | p2 | p3)
  dev.off()

  # rescue unassigned clusters: assign them to lineages based on lower resolution clustering
  new_marker_assignment = marker_assignment
  for (i in c(1:length(marker_assignment))){
    if (is.na(marker_assignment[i])){
      original_clusters = table(combined_seurat_subcluster$wsnn_res.4[combined_seurat_subcluster$wsnn_res.5.5 == i-1])
      original_clusters = names(original_clusters[original_clusters > sum(original_clusters) * 0.1])
      original_clusters_cell_types = table(combined_seurat_subcluster$cell_assignment[combined_seurat_subcluster$wsnn_res.4 %in% original_clusters])
      original_clusters_cell_types = names(original_clusters_cell_types[original_clusters_cell_types > sum(original_clusters_cell_types) * 0.1])
      original_clusters_cell_types = original_clusters_cell_types[!is.na(original_clusters_cell_types)]
      if (is.null(original_clusters_cell_types) || length(original_clusters_cell_types) == 0){
        next
      }
      original_clusters_cell_types = unlist(strsplit(original_clusters_cell_types, split = "early_"))[unlist(strsplit(original_clusters_cell_types, split = "early_")) != ""]
      original_clusters_cell_types = unique(str_replace_all(original_clusters_cell_types, "[a-z]", "x"))
      if (length(original_clusters_cell_types) == 1){
        new_marker_assignment[i] = original_clusters_cell_types
      }
    }
  }
  new_metadata = as.numeric(combined_seurat_subcluster$seurat_clusters)
  for (cluster_n in c(1:length(new_marker_assignment))){
    if (!is.na(new_marker_assignment[cluster_n])){
      new_metadata[as.numeric(combined_seurat_subcluster$seurat_clusters) == cluster_n] = new_marker_assignment[cluster_n]
    }
    else{
      new_metadata[as.numeric(combined_seurat_subcluster$seurat_clusters) == cluster_n] = NA
    }
  }
  combined_seurat_subcluster$cell_assignment = new_metadata
  p1 = DimPlot(combined_seurat_subcluster, reduction = "clusterumapWNN", order = T, label = T, group.by = "cell_assignment")
  p2 = DimPlot(combined_seurat_subcluster, reduction = "clusterumapWNN", order = T, label = T, group.by = "seurat_clusters") & NoLegend()
  p3 = DimPlot(combined_seurat_subcluster, reduction = "clusterumapWNN", order = T, label = T, group.by = "wsnn_res.4") & NoLegend()
  pdf(file = paste0(plot_dir_w_date, "/", lineage, "_full_annotation.UMAP.", annotation, ".pdf"), width = 14, height = 6, useDingbats = F)
  print(p1 | p2 | p3)
  dev.off()

  output_table = data.frame("cell_type" = combined_seurat_subcluster$cell_assignment, row.names = colnames(combined_seurat_subcluster))
  write.table(output_table, file = paste0(output_dir, "/", lineage, ".cell_assignment.", annotation, ".txt"), quote = F, sep = "\t", row.names = T, col.names = T)
}
