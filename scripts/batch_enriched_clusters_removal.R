# QC of the merged dataset and removal of batch-enriched clusters (adapted from batch_enriched_clusters_removal.20240430.R)
# changes:
#  * samples, batches and stage groups come from samples.tsv instead of 7 hard-coded experiments / 3 batches
#  * batch rule generalised: within a stage group, a cluster is flagged if one batch fraction > BATCH_FOLD x another
#    batch fraction (paper: early batches 1 vs 2, 4x); a cluster is also flagged if batches of LATE_STAGE_GROUPS
#    contribute > LATE_MAX_FRACTION of its cells (paper: batch 3 > 0.8)
#  * SoupX soupRange per batch = median barhop / ambient of its samples (paper: 100-400, 20-150, 100-200)
#  * total fragments matched by barcode; QC plotting code factorised
# L3 adaptation:
#  * QC uses an L3 germline marker score (tissueMarkers consensus, germline panel) instead of the embryonic
#    P-granule transcripts, and adds % ribosomal-protein UMIs; rRNA is reported in % like mito
#  * the embryo-specific protection of 2-cell / ABx clusters (pes-10, vet-2, cpg-3, ccch-2) is removed;
#    clusters can optionally be protected from removal with PROTECT_CLUSTER_MARKERS (config.sh)
#  * batch flagging corrected for batch size, groups only compared when they contribute meaningfully to a
#    cluster (BATCH_MIN_GROUP_NUCLEI / BATCH_MIN_GROUP_FRACTION), report-only mode (BATCH_REMOVAL_MODE),
#    per-cluster table batch_enrichment_table.tsv with the reason for flagging and top marker genes
#  * appends its filtering steps to the per-sample QC table qc_summary/<sample>.<annotation>.qc_steps.tsv
library(Seurat)
library(Signac)
library(ggplot2)
library(Hmisc)
library(SoupX)
library(glmGamPoi)
library(presto)
source("scripts/pipeline_utils.R")
# --- Seurat v4/v5 compatibility: manual module score using LayerData ---
# Replaces AddModuleScore which crashes on Assay5 (SeuratObject 5.0.0) via
# nested withCallingHandlers accessing @counts slot that doesn't exist in Assay5.
manual_module_score <- function(obj, genes, assay_name = "RNA", layer_name = "data", col_name = "score1") {
  expr <- LayerData(obj, assay = assay_name, layer = layer_name)
  genes_use <- intersect(genes, rownames(expr))
  message("module score '", col_name, "': using ", length(genes_use), " of ", length(genes), " features")
  if (length(genes_use) == 0) {
    obj[[col_name]] <- 0
    return(obj)
  }
  obj[[col_name]] <- as.numeric(Matrix::colMeans(expr[genes_use, , drop = FALSE]))
  obj
}


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
combined_seurat_UMAP = (DimPlot(combined_seurat, reduction = "umapWNN", label = T) +
  labs(title = "WNN UMAP – all clusters")) & NoLegend()

## optionally protect clusters marked by given genes from removal. The paper protected the 2-cell / ABx
## embryo clusters (cpg-3 / ccch-2), which were under-represented in some batches because of staging
## differences. Not needed for L3, so empty by default; set PROTECT_CLUSTER_MARKERS in config.sh if a rare
## population should never be removed (a cluster is protected if any listed gene is one of its positive markers).
protect_genes = strsplit(env_chr("PROTECT_CLUSTER_MARKERS", ""), "[ ,]+")[[1]]
protect_genes = intersect(protect_genes[nzchar(protect_genes)], rownames(combined_seurat[["RNA"]]))
protected_clusters = character(0)
if (length(protect_genes) > 0) {
  DefaultAssay(combined_seurat) = "RNA"
  protect_markers = tryCatch(
    FindAllMarkers(combined_seurat, assay = "RNA", features = protect_genes, logfc.threshold = 0.5, only.pos = TRUE),
    error = function(e) { warning("FindAllMarkers for protected genes failed: ", conditionMessage(e)); NULL })
  if (!is.null(protect_markers) && nrow(protect_markers) > 0)
    protected_clusters = unique(as.character(protect_markers$cluster))
  message("clusters protected from removal (", paste(protect_genes, collapse = ", "), "): ",
          if (length(protected_clusters)) paste(protected_clusters, collapse = ",") else "none")
}

## proportion of cells from each experiment / batch per cluster
cluster_levels = levels(combined_seurat$seurat_clusters)
exp_cell_mtx = t(sapply(exp_ids, function(id) table(factor(combined_seurat$seurat_clusters[combined_seurat$exp == id], levels = cluster_levels))))
rownames(exp_cell_mtx) = exp_ids
colnames(exp_cell_mtx) = cluster_levels
exp_cell_mtx_fraction = apply(exp_cell_mtx, 2, function(x){x/sum(x)})

batch_cell_mtx = rowsum(exp_cell_mtx, batch_of_sample[rownames(exp_cell_mtx)])
# share of each cluster's nuclei coming from each batch (as in the paper)
batch_cell_mtx_fraction = apply(batch_cell_mtx, 2, function(x){x/sum(x)})
if (is.null(dim(batch_cell_mtx_fraction))) batch_cell_mtx_fraction = matrix(batch_cell_mtx_fraction, nrow = 1, dimnames = list(rownames(batch_cell_mtx), cluster_levels))
# share of each batch's own nuclei that fall in each cluster: corrects for batches of different size
batch_totals = rowSums(batch_cell_mtx)
batch_prop = sweep(batch_cell_mtx, 1, batch_totals, "/")
# cluster composition as if every batch had the same number of nuclei (columns sum to 1)
batch_prop_norm = apply(batch_prop, 2, function(x) if (sum(x) > 0) x / sum(x) else x)
if (is.null(dim(batch_prop_norm))) batch_prop_norm = matrix(batch_prop_norm, nrow = 1, dimnames = list(rownames(batch_cell_mtx), cluster_levels))
cluster_n = colSums(batch_cell_mtx)

## flag clusters whose batch composition is skewed
# Within each stage group, the batches are compared on the share of their own nuclei in the cluster; the
# cluster is flagged if the highest share is > BATCH_FOLD x the lowest. A group is only compared for a
# cluster when it contributes meaningfully to it (>= BATCH_MIN_GROUP_NUCLEI nuclei AND >= BATCH_MIN_GROUP_FRACTION
# of the cluster), so a cluster specific to one group (e.g. one strain) is not flagged because of a few stray
# nuclei of the other group landing in one of its replicates.
min_group_fraction = env_num("BATCH_MIN_GROUP_FRACTION", 0.2)
min_group_nuclei = env_num("BATCH_MIN_GROUP_NUCLEI", 30)
removal_mode = env_chr("BATCH_REMOVAL_MODE", "remove")
if (removal_mode %nin% c("remove", "report")) stop("BATCH_REMOVAL_MODE must be 'remove' or 'report'")
flag_reason = setNames(rep("", length(cluster_levels)), cluster_levels)
add_reason = function(cl, txt) flag_reason[cl] <<- if (flag_reason[cl] == "") txt else paste(flag_reason[cl], txt, sep = "; ")
for (g in unique(batch_group)) {
  b = names(batch_group)[batch_group == g]
  if (length(b) < 2) next
  g_n = colSums(batch_cell_mtx[b, , drop = FALSE])
  present = g_n >= min_group_nuclei & g_n / pmax(cluster_n, 1) >= min_group_fraction
  for (cl in cluster_levels[present]) {
    p = batch_prop[b, cl]
    hi = which.max(p); lo = which.min(p)
    if (p[hi] > batch_fold * p[lo])
      add_reason(cl, paste0(g, ": ", b[hi], " vs ", b[lo], " = ",
                            if (p[lo] == 0) "only in one batch" else paste0(round(p[hi] / p[lo], 1), "x")))
  }
}
late_batches = names(batch_group)[batch_group %in% late_groups]
if (length(late_batches) > 0) {
  late_frac = colSums(batch_prop_norm[late_batches, , drop = FALSE])
  for (cl in cluster_levels[late_frac > late_max_fraction])
    add_reason(cl, paste0("stage groups ", paste(late_groups, collapse = ","), " make up ", round(100 * late_frac[cl]), "%"))
}
flagged_clusters = cluster_levels[flag_reason != ""]
# never remove protected clusters (PROTECT_CLUSTER_MARKERS)
batch_enriched_clusters = flagged_clusters[as.character(flagged_clusters) %nin% protected_clusters]
batch_removed_clusters = if (removal_mode == "remove") batch_enriched_clusters else character(0)
message("batch-enriched clusters: ", if (length(batch_enriched_clusters)) paste(batch_enriched_clusters, collapse = ",") else "none",
        " (BATCH_REMOVAL_MODE = ", removal_mode, if (removal_mode == "report") ": flagged clusters are kept" else "", ")")

combined_seurat$batch_enrichment = ifelse(combined_seurat$seurat_clusters %in% batch_removed_clusters, "flagged - removed",
                                   ifelse(combined_seurat$seurat_clusters %in% batch_enriched_clusters, "flagged - kept (report mode)",
                                   ifelse(combined_seurat$seurat_clusters %in% flagged_clusters, "flagged - protected", "not flagged")))
p_batch = DimPlot(combined_seurat, group.by = "batch_enrichment", order=T, reduction = "umapWNN") +
  labs(title = "Batch-enriched clusters",
       subtitle = paste0("BATCH_REMOVAL_MODE = ", removal_mode, "; details in batch_enrichment_table.tsv"))

pdf(file=paste0(plot_dir_w_date, "/merged_seurat.", annotation, ".batch_enrichment.pdf"), width=16, height = 8)
par(mfrow=c(1,1), mar = c(7, 4, 4, 2))
barplot(exp_cell_mtx_fraction, beside = F, las = 1, legend.text = rownames(exp_cell_mtx_fraction),
        main = "Share of each cluster's nuclei from each experiment (sample)",
        ylab = "Fraction of the cluster's nuclei", xlab = "Cluster")
barplot(batch_prop_norm, beside = F, las = 1, legend.text = rownames(batch_prop_norm),
        main = "Batch composition of each cluster, corrected for batch size (used for flagging)",
        ylab = "Fraction (each batch scaled to equal size)", xlab = "Cluster")
mtext(paste0("Within a stage group, a cluster is flagged if one batch's share of its own nuclei is > ", batch_fold,
             "x another's; only groups with >= ", min_group_nuclei, " nuclei and >= ", 100 * min_group_fraction,
             "% of the cluster are compared"), side = 1, line = 5, cex = 0.75)
if (length(late_batches) > 0) abline(h = late_max_fraction, lwd = 2, lty = 2)
print(combined_seurat_UMAP | p_batch)
dev.off()

DefaultAssay(combined_seurat) = "RNA"
combined_seurat = NormalizeData(combined_seurat)

## table of the batch composition of every cluster, with top marker genes of the flagged clusters
top_markers = setNames(rep("", length(cluster_levels)), cluster_levels)
Idents(combined_seurat) = "seurat_clusters"
for (cl in flagged_clusters) {
  m = tryCatch(FindMarkers(combined_seurat, ident.1 = cl, assay = "RNA", only.pos = TRUE, logfc.threshold = 0.5,
                           max.cells.per.ident = 500, verbose = FALSE),
               error = function(e) { warning("FindMarkers failed for cluster ", cl, ": ", conditionMessage(e)); NULL })
  if (!is.null(m) && nrow(m) > 0) {
    m = m[m$p_val_adj < 0.05, , drop = FALSE]
    top_markers[cl] = paste(head(rownames(m)[order(-m$avg_log2FC)], 10), collapse = ",")
  }
}
enrich_tab = data.frame(cluster = cluster_levels, n_nuclei = cluster_n, stringsAsFactors = FALSE)
for (b in rownames(batch_cell_mtx)) {
  enrich_tab[[paste0("n_", b)]] = batch_cell_mtx[b, ]
  enrich_tab[[paste0("pct_of_", b, "_nuclei")]] = round(100 * batch_prop[b, ], 2)
}
enrich_tab$flagged = cluster_levels %in% flagged_clusters
enrich_tab$protected = cluster_levels %in% flagged_clusters & cluster_levels %in% protected_clusters
enrich_tab$removed = cluster_levels %in% batch_removed_clusters
enrich_tab$reason = unname(flag_reason)
enrich_tab$top_markers_flagged = unname(top_markers)
write.table(enrich_tab, file = paste0(plot_dir_w_date, "/merged_seurat.", annotation, ".batch_enrichment_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# QC 1: L3 germline marker expression (paper: embryonic P-granule transcripts)
# germline panel of the L3 consensus marker table; labels translated to the gene names of this
# annotation via WBGeneID (some labels come from newer WormBase releases)
l3_markers = read.csv(env_chr("L3_MARKER_RANKED_FILE", "tissueMarkers/consensus/markers_consensus_min2sources.csv"),
                      stringsAsFactors = FALSE)
gene_id_name = read.table(paste0(canon, ".gene_name.txt"), col.names = c("wbgene", "gene_name"))
ws_name = setNames(gene_id_name$gene_name, gene_id_name$wbgene)[l3_markers$WBGeneID]
l3_markers$gene = ifelse(is.na(ws_name), l3_markers$gene_label, ws_name)
germline_genes = unique(l3_markers$gene[l3_markers$canonical_tissue == env_chr("GERMLINE_TISSUE", "germline") &
                                        !as.logical(l3_markers$multi_tissue_marker)])
# AddModuleScore replaced: crashes on Assay5 via withCallingHandlers @counts access
combined_seurat = manual_module_score(combined_seurat, germline_genes, col_name = "germline_score")

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
  # @counts on Assay5 replaced with rownames(combined_seurat[["RNA"]]) which is Assay5-safe
  raw_b = raw_b[row.names(raw_b) %in% rownames(combined_seurat[["RNA"]]),]
  # LayerData is the SeuratObject v5 API, safe on Assay5
  counts_b = LayerData(combined_seurat[, combined_seurat$exp %in% ids], assay = "RNA", layer = "counts")
  counts_b = counts_b[row.names(counts_b) %in% row.names(raw_b),]
  soup_range = c(median(as.numeric(samples$barhop[samples$batch == b])), median(as.numeric(samples$ambient[samples$batch == b])))
  sc_b = SoupChannel(raw_b, counts_b, calcSoupProfile=F)
  sc_b = estimateSoup(sc_b, soupRange=soup_range)
  row.names(sc_b$soupProfile[sc_b$soupProfile$est > 0.0005,])
})
high_ambient = Reduce(intersect, high_ambient_per_batch)
# AddModuleScore replaced: crashes on Assay5
combined_seurat = manual_module_score(combined_seurat, high_ambient, col_name = "ambient1")

# QC 3: rRNA, lincRNA and mito RNA expression
mt_genes <- read.table(paste0(canon, ".filtered.MtDNA_genes.genes"))
gene_type = read.table(paste0(canon, ".gene_name_gene_type.txt"))
rrna_genes <-  gene_type$V1[gene_type$V2  == "rRNA" & gene_type$V1 %nin% mt_genes$V1]
ribo_prot_genes = grep("^rp[ls]-[0-9]", rownames(raw_matrices[[1]]), value = TRUE)
# % of raw UMIs per barcode (same definition as seurat_soupx.R)
pct_in = function(genes) {
  v = unlist(lapply(raw_matrices, function(m) 100 * Matrix::colSums(m[row.names(m) %in% genes, , drop = FALSE]) /
                                              pmax(Matrix::colSums(m), 1)), use.names = FALSE)
  names(v) = unlist(lapply(raw_matrices, colnames), use.names = FALSE)
  v
}
combined_seurat$rrna_pct      = pct_in(rrna_genes)[colnames(combined_seurat)]
combined_seurat$ribo_prot_pct = pct_in(ribo_prot_genes)[colnames(combined_seurat)]
rm(raw_matrices); gc()

lincRNA_genes <- gene_type$V1[gene_type$V2  == "lincRNA" & gene_type$V1 %nin% mt_genes$V1]
# AddModuleScore replaced: crashes on Assay5
combined_seurat = manual_module_score(combined_seurat, lincRNA_genes, col_name = "lincRNA_RNA1")

# QC 4: FRiP / ATAC feature number
combined_seurat$tot_frag = total_fragments(combined_seurat, "ATAC")
combined_seurat = FRiP(object = combined_seurat, assay = "ATAC", total.fragments = "tot_frag", col.name = "FRiP")

qc_colours = c("grey90", "cyan", "deepskyblue", "darkorange", "darkred", "darkorchid3")

# Human-readable labels for each QC metric plotted by qc_plots()
QC_LABELS <- c(
  germline_score = "L3 germline marker score (germline nuclei / germline RNA in the soup)",
  ambient1     = "Ambient RNA score (SoupX background contamination)",
  mito_pct     = "% mitochondrial UMIs (cytoplasmic: ambient RNA / damaged nuclei)",
  rrna_pct     = "% rRNA UMIs (cytoplasmic: ambient RNA)",
  ribo_prot_pct = "% ribosomal-protein UMIs (cytoplasmic: ambient RNA)",
  lincRNA_RNA1 = "lincRNA module score",
  FRiP         = "Fraction of Reads in Peaks (ATAC quality)",
  nFeature_WNN = "WNN feature count (combined RNA + ATAC complexity)"
)
qc_label <- function(f) if (f %in% names(QC_LABELS)) QC_LABELS[[f]] else f

# qc_plots rewritten with ggplot2 to avoid FeaturePlot/VlnPlot touching @counts on Assay5
# (metadata features like germline_score are safe to plot but VlnPlot still traverses assay internals)
qc_plots = function(feature, log = FALSE, max.cutoff = NA) {
  umap_df <- as.data.frame(Embeddings(combined_seurat, reduction = "umapWNN"))
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  vals <- combined_seurat@meta.data[[feature]]
  if (!is.null(max.cutoff) && !is.na(max.cutoff)) vals <- pmin(vals, max.cutoff, na.rm = TRUE)
  df_umap <- cbind(umap_df, Value = vals)
  df_umap <- df_umap[order(df_umap$Value, na.last = FALSE), ]
  p_feat <- ggplot(df_umap, aes(x = UMAP_1, y = UMAP_2, colour = Value)) +
    geom_point(size = 0.3) +
    scale_colour_gradientn(colours = qc_colours, na.value = "grey90") +
    labs(title = qc_label(feature), subtitle = paste0("(", feature, ")")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8))
  print(combined_seurat_UMAP | p_feat)
  means <- tapply(combined_seurat@meta.data[[feature]], combined_seurat$seurat_clusters, mean, na.rm = TRUE)
  ordered_clusters <- names(sort(means, decreasing = TRUE))
  cols <- ifelse(ordered_clusters %in% batch_enriched_clusters, "darkred", "deepskyblue")
  df_vln <- data.frame(
    cluster = factor(combined_seurat$seurat_clusters, levels = ordered_clusters),
    value   = combined_seurat@meta.data[[feature]]
  )
  p_vln <- ggplot(df_vln, aes(x = cluster, y = value, fill = cluster)) +
    geom_violin(scale = "width") +
    scale_fill_manual(values = setNames(cols, ordered_clusters)) +
    labs(title = qc_label(feature),
         subtitle = "Clusters sorted by mean value  ·  darkred = batch-enriched (flagged for removal)  ·  deepskyblue = retained") +
    xlab("Cluster (sorted by mean value)") + ylab(feature) +
    theme_classic() + NoLegend() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8))
  if (log) p_vln <- p_vln + scale_y_log10()
  print(p_vln)
}

pdf(file=paste0(plot_dir_w_date, "/merged_seurat.", annotation, ".QC.pdf"), width=16, height = 8)
qc_plots("germline_score")
qc_plots("ambient1")
qc_plots("mito_pct")
qc_plots("rrna_pct")
qc_plots("ribo_prot_pct")
qc_plots("lincRNA_RNA1")
qc_plots("FRiP")
qc_plots("nFeature_WNN", log = TRUE, max.cutoff = 5000)
dev.off()

# remove batch-enriched cells and cells with few WNN features, then re-process each experiment and re-cluster
cells_to_keep = colnames(combined_seurat)[combined_seurat$nFeature_WNN > min_wnn_features & combined_seurat$seurat_clusters %nin% batch_removed_clusters]
message(length(cells_to_keep), " of ", ncol(combined_seurat), " nuclei retained")

# step-by-step QC summary: append the stage 7 rows to qc_summary/<sample>.<annotation>.qc_steps.tsv
md = combined_seurat@meta.data
atac_col = intersect(c("nCount_ATAC", "nCount_WNN"), colnames(md))[1]
pass_wnn = md$nFeature_WNN > min_wnn_features
pass_batch = pass_wnn & md$seurat_clusters %nin% batch_removed_clusters
none = function(x) if (length(x)) paste(x, collapse = ",") else "none"
for (id in exp_ids) {
  in_sample = md$exp == id
  r7 = function(step, keep, params, note = "")
    qc_row(7, step, params, md$nCount_RNA[keep], md$nFeature_RNA[keep],
           if (!is.na(atac_col)) md[[atac_col]][keep] else NULL,
           paste0(atac_col, " (ATAC counts in the bulk peak set)"), note)
  qc_write(list(
    r7("Merged object (stage 6)", in_sample, "input to stage 7",
       "RNA metrics from the SoupX-corrected counts"),
    r7("Low WNN feature filter", in_sample & pass_wnn,
       paste0("nFeature_WNN > MIN_WNN_FEATURES = ", min_wnn_features, " (peaks of the cluster-specific peak set with >= 1 count)")),
    r7("Batch-enriched cluster removal", in_sample & pass_batch,
       paste0("BATCH_REMOVAL_MODE = ", removal_mode, "; WNN clusters (resolution 5) flagged: ", none(batch_enriched_clusters),
              "; removed: ", none(batch_removed_clusters), "; BATCH_FOLD = ", batch_fold,
              " on batch-size-corrected shares; groups compared if >= ", min_group_nuclei, " nuclei and >= ",
              100 * min_group_fraction, "% of the cluster",
              "; LATE_MAX_FRACTION = ", late_max_fraction, " for stage groups ", paste(late_groups, collapse = ","),
              "; protected clusters: ", none(protected_clusters)),
       "needs >= 2 batches per stage group, or a stage group in LATE_STAGE_GROUPS, to flag anything")),
    qc_summary_file(id, annotation), replace_stage = "7")
}
# Use Signac Fragments() API instead of @fragments slot access
all_fragments = Fragments(combined_seurat[["WNN"]])

clean_objects = lapply(exp_ids, function(id) {
  keep = combined_seurat$exp == id & colnames(combined_seurat) %in% cells_to_keep
  # LayerData is SeuratObject v5 API, safe on Assay5 (avoids @counts slot)
  rna_mtx = LayerData(combined_seurat[, keep], assay = "RNA", layer = "counts")
  # WNN is ChromatinAssay (Signac); use LayerData for consistency
  atac_mtx = LayerData(combined_seurat[, keep], assay = "WNN", layer = "counts")
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

# @scale.data slot does not exist in Assay5; use LayerData (SeuratObject v5 API)
VariableFeatures(gcdata[["SCT"]]) <- rownames(LayerData(gcdata, assay = "SCT", layer = "scale.data"))
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

# PrepSCTFindMarkers() removed: FindMarkers is never called in this script,
# and it fails on merged Assay5 (SCT v2) objects in Seurat v5.

# Join split per-sample Assay5 layers (Seurat v5 / SeuratObject 5.0.0 compatibility)
# Without this, FeaturePlot crashes with "no slot of name counts for Assay5"
# JoinLayers only on plain Assay5 (RNA) — SCTAssay and ChromatinAssay lack this method
if ("RNA" %in% Assays(gcdata) && inherits(gcdata[["RNA"]], "Assay5")) {
  gcdata[["RNA"]] <- JoinLayers(gcdata[["RNA"]])
}

p1 <- DimPlot(gcdata, reduction = "umap", group.by = "wsnn_res.4", label = TRUE, label.size = 4) +
      labs(title = "RNA-only UMAP", subtitle = "SCT normalisation, 50 PCA dims")
p2 <- DimPlot(gcdata, reduction = "umapATAC", group.by = "wsnn_res.4", label = TRUE, label.size = 4) +
      labs(title = "ATAC-only UMAP", subtitle = "LSI dimensionality reduction, dims 2-50")
p3 <- DimPlot(gcdata, reduction = "umapWNN", group.by = "wsnn_res.4", label = TRUE, label.size = 4) +
      labs(title = "Weighted Nearest Neighbour UMAP", subtitle = "Integrated RNA + ATAC (WNN)")
pdf(file=paste0(plot_dir_w_date, "/all_samples_SCT_clean.", annotation, ".umap_WNN.pdf"), width=24, height = 8)
print(p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5),
                                          plot.subtitle = element_text(hjust = 0.5, size = 9)))
dev.off()

pdf(file=paste0(plot_dir_w_date, "/all_samples_SCT_clean.", annotation, ".umap_split.WNN.pdf"), width=24, height = 8)
p1 = DimPlot(gcdata, reduction = "umapWNN", label = TRUE, label.size = 4) +
     labs(title = "WNN UMAP – cluster labels")
p2 = DimPlot(gcdata, reduction = "umapWNN", group.by = "exp") +
     labs(title = "WNN UMAP – coloured by sample")
print(p1 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)) | p2  & theme(plot.title = element_text(hjust = 0.5)))
print(DimPlot(gcdata, reduction = "umapWNN", split.by = "exp") & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
dev.off()

# L3 key tissue marker plots (top 3 per broad tissue from markers_shared_min2studies.tsv,
# i.e. markers supported by >= 2 published C. elegans studies).
# Uses LayerData() + Embeddings() + ggplot2 directly, avoiding FeaturePlot's
# internal withCallingHandlers which crashes on merged Assay5 (Seurat v4/v5 mismatch).
# L3 key tissue markers - data-driven selection across 6 studies
# Ranked by: (1) #studies with strict specificity (other_max==0, tissue_max>=1)
#            (2) total support, (3) preferring named genes
l3_key_markers <- c(
  "inos-1",  "aman-1",  "clec-145",   # Coelomocytes  (2 strict studies)
  "cpr-3",   "klo-1",   "nspc-1",     # Excretory     (1 strict study)
  "pgl-1",   "pid-1",   "sip-1",      # Germline      (4 strict studies)
  "glt-1",   "grl-12",  "cka-1",      # Glia          (2 strict studies)
  "inx-8",   "inx-9",   "C47D12.5",   # Gonad         (2 strict studies)
  "col-107", "dpy-13",  "ram-2",      # Hypodermis    (2 strict studies)
  "asp-1",   "act-5",   "cpr-9",      # Intestine     (5 strict studies; F57F5.1 = cpr-9 in WS285)
  "ace-1",   "acs-16",  "bgal-1",     # Muscle        (1 strict study)
  "acc-1",   "acr-2",   "acr-17",     # Neurons       (2 strict studies)
  "marg-1",  "abf-6",   "abu-14"      # Pharynx       (2 strict studies)
)
l3_marker_tissue <- c(
  "inos-1"="Coelomocytes",  "aman-1"="Coelomocytes",  "clec-145"="Coelomocytes",
  "cpr-3"="Excretory",      "klo-1"="Excretory",      "nspc-1"="Excretory",
  "pgl-1"="Germline",       "pid-1"="Germline",        "sip-1"="Germline",
  "glt-1"="Glia",           "grl-12"="Glia",           "cka-1"="Glia",
  "inx-8"="Gonad",          "inx-9"="Gonad",           "C47D12.5"="Gonad",
  "col-107"="Hypodermis",   "dpy-13"="Hypodermis",     "ram-2"="Hypodermis",
  "asp-1"="Intestine",      "act-5"="Intestine",       "cpr-9"="Intestine",
  "ace-1"="Muscle",         "acs-16"="Muscle",         "bgal-1"="Muscle",
  "acc-1"="Neurons",        "acr-2"="Neurons",         "acr-17"="Neurons",
  "marg-1"="Pharynx",       "abf-6"="Pharynx",         "abu-14"="Pharynx"
)

message("Plotting L3 key tissue markers using LayerData + ggplot2 (Seurat v4/v5 compat)")
umap_df <- as.data.frame(Embeddings(gcdata, reduction = "umapWNN"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")

# Use RNA normalized data (NormalizeData was called earlier); LayerData is SeuratObject v5 API
rna_layer <- LayerData(gcdata, assay = "RNA", layer = "data")
genes_present <- intersect(l3_key_markers, rownames(rna_layer))
message("L3 key markers found in object: ", length(genes_present), " / ", length(l3_key_markers))

pdf(file = paste0(plot_dir_w_date, "/all_samples_SCT.", annotation, ".L3_key_markers.WNN_clean.pdf"),
    width = 8, height = 8)
for (g in genes_present) {
  expr <- as.numeric(rna_layer[g, ])
  df   <- cbind(umap_df, Expression = expr)
  df   <- df[order(df$Expression), ]   # low expr plotted first so high expr on top
  tissue_label <- l3_marker_tissue[g]
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, colour = Expression)) +
    geom_point(size = 0.4) +
    scale_colour_gradientn(colours = feature_colours) +
    ggtitle(paste0(g, "  [", tissue_label, "]")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
}
dev.off()
message("L3 key marker plots done")

saveRDS(gcdata, file = out_rds)
write.table(gcdata[[]], file = sub("\\.rds$", ".table.txt", out_rds), quote = F, sep = "\t", row.names = T, col.names = T)
