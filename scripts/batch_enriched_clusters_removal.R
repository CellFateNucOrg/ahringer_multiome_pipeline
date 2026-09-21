# QC of the merged dataset and removal of batch-enriched clusters (adapted from batch_enriched_clusters_removal.20240430.R)
# changes:
#  * samples, batches and stage groups come from samples.tsv instead of 7 hard-coded experiments / 3 batches
#  * batch rule generalised: within a stage group, a cluster is flagged if one batch fraction > BATCH_FOLD x another
#    batch fraction (paper: early batches 1 vs 2, 4x); a cluster is also flagged if batches of LATE_STAGE_GROUPS
#    contribute > LATE_MAX_FRACTION of its cells (paper: batch 3 > 0.8)
#  * SoupX soupRange per batch = median barhop / ambient of its samples (paper: 100-400, 20-150, 100-200)
#  * total fragments matched by barcode; QC plotting code factorised
library(Seurat)
library(Signac)
library(ggplot2)
library(Hmisc)
library(SoupX)
library(glmGamPoi)
library(presto)
source("scripts/pipeline_utils.R")

args = commandArgs(trailingOnly=TRUE)
plot_dir_w_date = make_plot_dir("batch_removal")

seurat_object = args[1]
annotation = args[2]
out_rds = sub("\\.WNN\\.rds$", ".WNN_clean.rds", seurat_object)
canon = canon_prefix()

samples = read_samples()
exp_ids = samples$sample_id
batch_of_sample = setNames(samples$batch, samples$sample_id)
batch_group = tapply(samples$stage_group, samples$batch, function(x) unique(x))
if (any(lengths(batch_group) != 1)) stop("each batch must belong to a single stage_group in samples.tsv")
batch_group = unlist(batch_group)
batch_fold = env_num("BATCH_FOLD", 4)
late_groups = strsplit(env_chr("LATE_STAGE_GROUPS", "late"), "[ ,]+")[[1]]
late_max_fraction = env_num("LATE_MAX_FRACTION", 0.8)
min_wnn_features = env_num("MIN_WNN_FEATURES", 200)

combined_seurat = readRDS(seurat_object)

### identify clusters showing a strong skew in their batch of origin
combined_seurat <- FindClusters(combined_seurat, graph.name = "wsnn", algorithm = 3, verbose = T, resolution = 5)
combined_seurat_UMAP = DimPlot(combined_seurat, reduction = "umapWNN", label = T) & NoLegend()

## mark the very early cells (2-cell and ABx) to avoid their removal (under-represented due to differences in staging)
DefaultAssay(combined_seurat) = "RNA"
all_markers = FindAllMarkers(combined_seurat,assay = "RNA", logfc.threshold = 0.5, only.pos = T)
early_clusters = all_markers[all_markers$gene %in% c("pes-10", "vet-2"),6]
twocell_cluster = all_markers[all_markers$gene == "cpg-3",6]
abx_cluster = all_markers[all_markers$gene == "ccch-2",6][which(all_markers[all_markers$gene == "ccch-2",6] %nin% early_clusters)]

## proportion of cells from each experiment / batch per cluster
cluster_levels = levels(combined_seurat$seurat_clusters)
exp_cell_mtx = t(sapply(exp_ids, function(id) table(factor(combined_seurat$seurat_clusters[combined_seurat$exp == id], levels = cluster_levels))))
rownames(exp_cell_mtx) = exp_ids
colnames(exp_cell_mtx) = cluster_levels
exp_cell_mtx_fraction = apply(exp_cell_mtx, 2, function(x){x/sum(x)})

batch_cell_mtx = rowsum(exp_cell_mtx, batch_of_sample[rownames(exp_cell_mtx)])
batch_cell_mtx_fraction = apply(batch_cell_mtx, 2, function(x){x/sum(x)})
if (is.null(dim(batch_cell_mtx_fraction))) batch_cell_mtx_fraction = matrix(batch_cell_mtx_fraction, nrow = 1, dimnames = list(rownames(batch_cell_mtx), cluster_levels))

flagged = rep(FALSE, ncol(batch_cell_mtx_fraction))
for (g in unique(batch_group)) {
  b = names(batch_group)[batch_group == g]
  if (length(b) < 2) next
  for (b1 in b) for (b2 in b) if (b1 != b2) {
    flagged = flagged | batch_cell_mtx_fraction[b1, ] > batch_cell_mtx_fraction[b2, ] * batch_fold
  }
}
late_batches = names(batch_group)[batch_group %in% late_groups]
if (length(late_batches) > 0) {
  flagged = flagged | colSums(batch_cell_mtx_fraction[late_batches, , drop = FALSE]) > late_max_fraction
}
batch_enriched_clusters = colnames(batch_cell_mtx_fraction)[flagged]
# never remove 2-cell / ABx clusters
batch_enriched_clusters = batch_enriched_clusters[as.character(batch_enriched_clusters) %nin% c(as.character(twocell_cluster), as.character(abx_cluster))]
message("batch-enriched clusters: ", paste(batch_enriched_clusters, collapse = ","))

combined_seurat$batch_enrichment = ifelse(combined_seurat$seurat_clusters %in% batch_enriched_clusters, "batch_enriched", "not_enriched")
p_batch = DimPlot(combined_seurat, group.by = "batch_enrichment", order=T, reduction = "umapWNN")

pdf(file=paste0(plot_dir_w_date, "/merged_seurat.", annotation, ".batch_enrichment.pdf"), width=16, height = 8)
par(mfrow=c(1,1))
barplot(exp_cell_mtx_fraction, beside = F, las = 1, main = "exp of origin")
abline(h=0.80, lwd = 2, lty = 2); abline(h=0.20, lwd = 2, lty = 2)
barplot(batch_cell_mtx_fraction, beside = F, las = 1, main = "batch of origin")
abline(h=0.80, lwd = 2, lty = 2); abline(h=0.20, lwd = 2, lty = 2)
print(combined_seurat_UMAP | p_batch)
dev.off()

DefaultAssay(combined_seurat) = "RNA"
combined_seurat = NormalizeData(combined_seurat)

# QC 1: P-granule gene expression
P_granule_genes = read.table("data/external_data/P_granule_transcripts.sorted.txt")
combined_seurat = AddModuleScore(combined_seurat, features = list(P_granule_genes$V1), name = "Pgranule")

# QC 2: ambient RNA - genes with high soup expression in every batch
raw_matrices = lapply(exp_ids, function(id) {
  m <- Seurat::Read10X(paste0("data/sc/rnaseq/", id, "/star_", annotation, "/star_", annotation, "_Solo.out/GeneFull/raw/um_reads"))
  colnames(m) <- paste0(id, "_", colnames(m))
  m
})
names(raw_matrices) = exp_ids

high_ambient_per_batch = lapply(unique(samples$batch), function(b) {
  ids = samples$sample_id[samples$batch == b]
  raw_b = do.call(cbind, raw_matrices[ids])
  raw_b = raw_b[row.names(raw_b) %in% row.names(combined_seurat@assays$RNA@counts),]
  counts_b = combined_seurat[,combined_seurat$exp %in% ids]@assays$RNA@counts
  counts_b = counts_b[row.names(counts_b) %in% row.names(raw_b),]
  soup_range = c(median(as.numeric(samples$barhop[samples$batch == b])), median(as.numeric(samples$ambient[samples$batch == b])))
  sc_b = SoupChannel(raw_b, counts_b, calcSoupProfile=F)
  sc_b = estimateSoup(sc_b, soupRange=soup_range)
  row.names(sc_b$soupProfile[sc_b$soupProfile$est > 0.0005,])
})
high_ambient = Reduce(intersect, high_ambient_per_batch)
combined_seurat = AddModuleScore(combined_seurat, features = list(high_ambient), name = "ambient")

# QC 3: rRNA, lincRNA and mito RNA expression
mt_genes <- read.table(paste0(canon, ".filtered.MtDNA_genes.genes"))
gene_type = read.table(paste0(canon, ".gene_name_gene_type.txt"))
rrna_genes <-  gene_type$V1[gene_type$V2  == "rRNA" & gene_type$V1 %nin% mt_genes$V1]
rrna_pct = unlist(lapply(raw_matrices, function(m) colSums(m[row.names(m) %in% rrna_genes, , drop = FALSE])/colSums(m)), use.names = FALSE)
names(rrna_pct) = unlist(lapply(raw_matrices, colnames), use.names = FALSE)
combined_seurat$rrna_pct = rrna_pct[colnames(combined_seurat)]
rm(raw_matrices); gc()

lincRNA_genes <- gene_type$V1[gene_type$V2  == "lincRNA" & gene_type$V1 %nin% mt_genes$V1]
combined_seurat = AddModuleScore(combined_seurat, features = list(c(lincRNA_genes)), assay = "RNA", name = "lincRNA_RNA", nbin = 10)

# QC 4: FRiP / ATAC feature number
combined_seurat$tot_frag = total_fragments(combined_seurat, "ATAC")
combined_seurat = FRiP(object = combined_seurat, assay = "ATAC", total.fragments = "tot_frag", col.name = "FRiP")

qc_colours = c("grey90", "cyan", "deepskyblue", "darkorange", "darkred", "darkorchid3")
qc_plots = function(feature, log = FALSE, max.cutoff = NA) {
  p = FeaturePlot(combined_seurat, label = T, repel = T, label.size = 2, label.color = 2, features = feature,
                  reduction = "umapWNN", order = T, max.cutoff = max.cutoff) & scale_colour_gradientn(colours = qc_colours, na.value = "grey90")
  print(combined_seurat_UMAP | p)
  means = tapply(combined_seurat@meta.data[[feature]], combined_seurat$seurat_clusters, mean, na.rm = TRUE)
  ordered_clusters = names(sort(means, decreasing = TRUE))
  cols = ifelse(ordered_clusters %in% batch_enriched_clusters, "darkred", "deepskyblue")
  print(VlnPlot(combined_seurat, features = feature, sort = "increasing", log = log, cols = cols) & NoLegend())
}

pdf(file=paste0(plot_dir_w_date, "/merged_seurat.", annotation, ".QC.pdf"), width=16, height = 8)
qc_plots("Pgranule1")
qc_plots("ambient1")
qc_plots("mito_pct")
qc_plots("rrna_pct")
qc_plots("lincRNA_RNA1")
qc_plots("FRiP")
qc_plots("nFeature_WNN", log = TRUE, max.cutoff = 5000)
dev.off()

# remove batch-enriched cells and cells with few WNN features, then re-process each experiment and re-cluster
cells_to_keep = colnames(combined_seurat)[combined_seurat$nFeature_WNN > min_wnn_features & combined_seurat$seurat_clusters %nin% batch_enriched_clusters]
message(length(cells_to_keep), " of ", ncol(combined_seurat), " nuclei retained")
all_fragments = combined_seurat@assays$WNN@fragments

clean_objects = lapply(exp_ids, function(id) {
  keep = combined_seurat$exp == id & colnames(combined_seurat) %in% cells_to_keep
  rna_mtx = combined_seurat[, keep]@assays$RNA@counts
  atac_mtx = combined_seurat[, keep]@assays$WNN@counts
  frag = all_fragments[[which(sapply(all_fragments, function(f) any(startsWith(Cells(f), paste0(id, "_")))))[1]]]
  o <- CreateSeuratObject(counts = rna_mtx, project = paste0(id, "_clean"))
  o <- SCTransform(o, vst.flavor = "v2", method = "glmGamPoi")
  o[["WNN"]] <- CreateChromatinAssay(counts = atac_mtx, sep = c(":", "-"), fragments = frag)
  DefaultAssay(o) = "WNN"
  o <- RunTFIDF(o)
  o@meta.data$exp <- id
  o
})

gcdata <- merge_seurat_list(clean_objects, merge.data = T)
rm(clean_objects); gc()

VariableFeatures(gcdata[["SCT"]]) <- rownames(gcdata[["SCT"]]@scale.data)
DefaultAssay(gcdata) = "SCT"
gcdata <- RunPCA(gcdata, seed.use = 42, features = VariableFeatures(object = gcdata), npcs = 100)
gcdata <- FindNeighbors(gcdata, dims = 1:50)
gcdata <- FindClusters(gcdata, resolution = 3, random.seed = 0)
gcdata <- RunUMAP(gcdata, seed.use = 42, dims = 1:50, return.model=TRUE)

DefaultAssay(gcdata) = "RNA"
gcdata = NormalizeData(gcdata)

DefaultAssay(gcdata) = "WNN"
gcdata <- FindTopFeatures(gcdata, min.cutoff = 'q0')
gcdata <- RunSVD(gcdata)
gcdata <- RunUMAP(gcdata, reduction = 'lsi', dims = 2:50, reduction.name = "umapATAC")
gcdata <- FindNeighbors(gcdata, reduction = 'lsi', dims = 2:50)
gcdata <- FindClusters(gcdata, verbose = T, algorithm = 3, resolution = 3)

gcdata <- FindMultiModalNeighbors(gcdata, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
gcdata <- RunUMAP(gcdata, nn.name = "weighted.nn", reduction.name = "umapWNN")
gcdata <- FindClusters(gcdata, graph.name = "wsnn", algorithm = 3, verbose = T, resolution = 4)

gcdata <- PrepSCTFindMarkers(gcdata)

p1 <- DimPlot(gcdata, reduction = "umap", group.by = "wsnn_res.4", label = TRUE, label.size = 4)
p2 <- DimPlot(gcdata, reduction = "umapATAC", group.by = "wsnn_res.4", label = TRUE, label.size = 4)
p3 <- DimPlot(gcdata, reduction = "umapWNN", group.by = "wsnn_res.4", label = TRUE, label.size = 4)
pdf(file=paste0(plot_dir_w_date, "/all_samples_SCT_clean.", annotation, ".umap_WNN.pdf"), width=24, height = 8)
print(p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
dev.off()

pdf(file=paste0(plot_dir_w_date, "/all_samples_SCT_clean.", annotation, ".umap_split.WNN.pdf"), width=24, height = 8)
p1 = DimPlot(gcdata, reduction = "umapWNN", label = TRUE, label.size = 4)
p2 = DimPlot(gcdata, reduction = "umapWNN", group.by = "exp")
print(p1 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)) | p2  & theme(plot.title = element_text(hjust = 0.5)))
print(DimPlot(gcdata, reduction = "umapWNN", split.by = "exp") & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
dev.off()

DefaultAssay(gcdata) = "SCT"
plot_key_markers(gcdata, "umapWNN", paste0(plot_dir_w_date, "/all_samples_SCT.", annotation, ".key_lineage_markers.WNN_clean.pdf"))

saveRDS(gcdata, file = out_rds)
write.table(gcdata[[]], file = sub("\\.rds$", ".table.txt", out_rds), quote = F, sep = "\t", row.names = T, col.names = T)
