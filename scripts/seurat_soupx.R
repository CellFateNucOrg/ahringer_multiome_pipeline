# Per-sample QC: mito / ATAC fragment filters, scDblFinder, SoupX, ATAC assay on the bulk peak set
# (adapted from seurat_soupx.20240207.R; commented-out exploratory code removed)
# changes: thresholds from samples.tsv / config, plots in plots/seurat_soupx/<sample>,
#          ATAC_count matched by barcode, fallback when P-cell binarisation is not possible
# L3 adaptation (replaces the embryo-specific gene sets of the paper):
#   * germline clusters excluded from the SoupX estimate are found with L3 germline markers from the
#     tissueMarkers consensus table (was: 492 embryonic P-granule transcripts)
#   * soup QC: mitochondrial, rRNA and ribosomal-protein UMI fractions compared between empty droplets
#     (soup) and nuclei. Mito transcripts are cytoplasmic, so in nuclei they mostly report soup; rRNA and
#     ribosomal-protein transcripts are also made in the nucleus, so they are shown for comparison only
#   * the SoupX soup profile is written to seurat_objects/<sample>.<annotation>.soupProfile.tsv, annotated
#     with the L3 tissue of each gene; the highest-ranking tissue-specific soup genes replace vet-6
#     (embryonic) in the decontamination change maps
#   * post-SoupX marker plots use the top L3 markers per tissue (was: embryonic lineage markers),
#     shown before and after SoupX
#   * step-by-step QC table qc_summary/<sample>.<annotation>.qc_steps.tsv: nuclei left after each cleaning step,
#     RNA UMIs / genes / ATAC counts per nucleus and the thresholds used (stage 7 appends its own steps)
library(ggplot2)
library(dplyr)
library(Seurat)
library(Signac)
library(patchwork)
library(DropletUtils)
library(SoupX)
library(pheatmap)
library(RColorBrewer)
library(scDblFinder)
library(scater)
library(Hmisc)
library(rtracklayer)
library(glmGamPoi)
library(multimode)

source("scripts/second_derivative_binarisation.R")
source("scripts/pipeline_utils.R")

args = commandArgs(trailingOnly=TRUE)

sample_id <- args[1]
annotation <- args[2]
min_frags = as.numeric(args[3])
merged_peaks_file = args[4]
samples = read_samples()
barhop_rna <- sample_value(samples, sample_id, "barhop")
lower_rna <- sample_value(samples, sample_id, "ambient")
mt_pct_threshold = env_num("MT_PCT_MAX", 5)
max_rna_umi = env_num("MAX_RNA_UMI", 10000)
soup_quantile = env_num("SOUPX_QUANTILE", 0.8)
canon = canon_prefix()

# L3 marker settings (config/config.sh)
l3_marker_file    = env_chr("L3_MARKER_RANKED_FILE", "tissueMarkers/consensus/markers_consensus_min2sources.csv")
germline_tissue   = env_chr("GERMLINE_TISSUE", "germline")
n_markers_plot    = env_num("N_MARKERS_PER_TISSUE_PLOT", 3)
n_soup_change     = env_num("N_SOUP_CHANGE_GENES", 6)
n_soup_score      = env_num("N_SOUP_SCORE_GENES", 50)

ED_barcodes = paste0("data/sc/rnaseq/", sample_id, "/ED_RNA_", annotation, "/ED_RNA_barcodes.txt")
set.seed(42)

plot_dir_w_date = make_plot_dir(file.path("seurat_soupx", sample_id))
seurat_objects_dir = "seurat_objects/"
dir.create(seurat_objects_dir, showWarnings = FALSE)

### L3 tissue markers ----------------------------------------------------------------------------
# ranked consensus table: canonical_tissue, consensus_rank, gene_label, WBGeneID, multi_tissue_marker, ...
# gene labels are translated to the WS285 names used in the count matrix via WBGeneID
# (some labels come from newer WormBase releases, e.g. spig-9, nssp-*)
l3_markers = read.csv(l3_marker_file, stringsAsFactors = FALSE)
gene_id_name = read.table(paste0(canon, ".gene_name.txt"), col.names = c("wbgene", "gene_name"))
ws_name = setNames(gene_id_name$gene_name, gene_id_name$wbgene)[l3_markers$WBGeneID]
l3_markers$gene = ifelse(is.na(ws_name), l3_markers$gene_label, ws_name)
l3_markers$multi_tissue_marker = as.logical(l3_markers$multi_tissue_marker)
l3_markers = l3_markers[order(l3_markers$canonical_tissue, l3_markers$consensus_rank), ]
# single-tissue markers only, for tissue look-ups and plots
l3_single = l3_markers[!l3_markers$multi_tissue_marker, ]
l3_single = l3_single[!(l3_single$gene %in% l3_single$gene[duplicated(l3_single$gene)]), ]
tissue_of_gene = setNames(l3_single$canonical_tissue, l3_single$gene)
# display names for plots / tables: the common name wherever one exists (current WormBase label if it is a
# common name, else the WS285 name, e.g. Y75B7AL.2 -> amn-1, F57F5.1 -> cpr-9); matching always uses WS285 names
l3_markers$display = ifelse(is_common_name(l3_markers$gene_label) | !is_common_name(l3_markers$gene),
                            l3_markers$gene_label, l3_markers$gene)
gene_display = with(l3_markers[!duplicated(l3_markers$gene), ], setNames(display, gene))
disp = function(g) unname(ifelse(g %in% names(gene_display), gene_display[g], g))
l3_tissues = unique(l3_markers$canonical_tissue)
pretty_tissue = function(t) gsub("_", " ", t)
germline_genes_all = l3_single$gene[l3_single$canonical_tissue == germline_tissue]
message("L3 markers: ", nrow(l3_single), " single-tissue genes in ", length(l3_tissues), " tissues; ",
        length(germline_genes_all), " ", germline_tissue, " genes")

raw_matrix_fran = paste0("data/sc/rnaseq/",sample_id, "/star_", annotation, "/star_", annotation, "_Solo.out/GeneFull/raw/um_reads")
rnaseq_data_raw <- Seurat::Read10X(raw_matrix_fran)

# read in the barcodes obtained by EmptyDrops and subset the object
combined_barcodes <- read.table(ED_barcodes)
rnaseq_data <- rnaseq_data_raw[,colnames(rnaseq_data_raw) %in% combined_barcodes$V1]
combined_data_seurat <- CreateSeuratObject(counts = rnaseq_data, project = sample_id )

### soup marker genes: mitochondrial, rRNA, ribosomal protein ----------------------------------------
mt_genes <- read.table(paste0(canon, ".filtered.MtDNA_genes.genes"))
gene_type = read.table(paste0(canon, ".gene_name_gene_type.txt"))
rrna_genes = gene_type$V1[gene_type$V2 == "rRNA" & gene_type$V1 %nin% mt_genes$V1]
ribo_prot_genes = grep("^rp[ls]-[0-9]", rownames(rnaseq_data_raw), value = TRUE)

raw_umi = Matrix::colSums(rnaseq_data_raw)
pct_of = function(genes) 100 * Matrix::colSums(rnaseq_data_raw[rownames(rnaseq_data_raw) %in% genes, , drop = FALSE]) / pmax(raw_umi, 1)
mt_pct_raw   = pct_of(mt_genes$V1)
rrna_pct_raw = pct_of(rrna_genes)
ribo_pct_raw = pct_of(ribo_prot_genes)
mt_genes_expression_pct = mt_pct_raw / 100   # fraction, as in the paper

# empty droplets in the SoupX soup range vs nuclei called by EmptyDrops
droplet_class = rep(NA_character_, length(raw_umi))
droplet_class[raw_umi >= as.numeric(barhop_rna) & raw_umi < as.numeric(lower_rna)] = "empty droplets = soup"
droplet_class[colnames(rnaseq_data_raw) %in% combined_barcodes$V1] = "nuclei (called by EmptyDrops)"
keep_idx = which(!is.na(droplet_class))
empty_idx = keep_idx[droplet_class[keep_idx] == "empty droplets = soup"]
if (length(empty_idx) > 20000) keep_idx = c(setdiff(keep_idx, empty_idx), sample(empty_idx, 20000))
soup_qc_df = data.frame(
  class  = rep(droplet_class[keep_idx], 3),
  metric = rep(c("% mitochondrial UMIs", "% rRNA UMIs", "% ribosomal protein UMIs"), each = length(keep_idx)),
  value  = c(mt_pct_raw[keep_idx], rrna_pct_raw[keep_idx], ribo_pct_raw[keep_idx]))
for (m in unique(soup_qc_df$metric)) for (cl in unique(soup_qc_df$class))
  message(sample_id, ": median ", m, " in ", cl, " = ", round(median(soup_qc_df$value[soup_qc_df$metric == m & soup_qc_df$class == cl]), 2))

pdf(file=paste0(plot_dir_w_date, "/", sample_id, ".", annotation, ".soup_marker_fractions.pdf"), width=12, height = 6)
print(ggplot(soup_qc_df, aes(x = class, y = value, fill = class)) +
      geom_violin(scale = "width") + geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white") +
      facet_wrap(~metric, scales = "free_y") +
      scale_fill_manual(values = c("empty droplets = soup" = "grey70", "nuclei (called by EmptyDrops)" = "deepskyblue")) +
      labs(title = paste0(sample_id, ": mitochondrial, ribosomal-protein and rRNA UMIs in empty droplets vs nuclei"),
           subtitle = paste0("Grey: empty droplets = barcodes with ", barhop_rna, "-", lower_rna,
                             " UMIs, i.e. the soup (same range SoupX uses). Blue: barcodes called as nuclei by EmptyDrops.\n",
                             "Mito transcripts are cytoplasmic, so in nuclei they mostly come from the soup and nuclei should have lower % than empty droplets.\n",
                             "rRNA and ribosomal-protein transcripts are also made in the nucleus (nucleolus / highly transcribed genes), so similar values are expected."),
           x = NULL, y = "% of UMIs in droplet") +
      theme_classic() + NoLegend() + theme(plot.subtitle = element_text(size = 8)))
dev.off()

### inspect mitochondrial contamination
sce.rnaseq_data <- as.SingleCellExperiment(combined_data_seurat)
is.mito <- rownames(sce.rnaseq_data) %in% mt_genes$V1
sce.rnaseq_data <- addPerCellQC(sce.rnaseq_data, subsets=list(Mt=is.mito))
print(summary(sce.rnaseq_data$subsets_Mt_percent == 0))

stats <- quickPerCellQC(colData(sce.rnaseq_data), sub.fields="subsets_Mt_percent")

pdf(file=paste0(plot_dir_w_date, "/", sample_id, ".", annotation, ".mitochondrial_pct_all_barcodes.pdf"), width=8, height = 8)
print(plotColData(sce.rnaseq_data, y="subsets_Mt_percent", colour_by=I(stats$high_subsets_Mt_percent)) +
      geom_hline(yintercept = mt_pct_threshold, linetype = 2) +
      labs(title = paste0(sample_id, ": % mitochondrial UMIs per nucleus (EmptyDrops barcodes)"),
           subtitle = paste0("Dashed line = MT_PCT_MAX (", mt_pct_threshold, "%); nuclei above it are removed. ",
                             "Colour = quickPerCellQC outlier call (informative only)")))
dev.off()

high_mt_cells = colnames(sce.rnaseq_data[,sce.rnaseq_data$subsets_Mt_percent >= mt_pct_threshold])

# remove high-Mito cells
rnaseq_data <- rnaseq_data[,colnames(rnaseq_data) %nin% high_mt_cells]

# remove MT genes, rRNAs and small RNAs from matrix
genes_to_retain = gene_type$V1[gene_type$V2 %in% c("lincRNA", "protein_coding", "pseudogene") & gene_type$V1 %nin% mt_genes$V1]
rnaseq_data <- rnaseq_data[row.names(rnaseq_data) %in% genes_to_retain,]

### step-by-step QC summary: barcodes left after each filter; RNA metrics on the genes kept for analysis
### (lincRNA, protein-coding, pseudogene; no MT genes) so that all steps are comparable
qc_mat = rnaseq_data_raw[row.names(rnaseq_data_raw) %in% genes_to_retain, ]
qc_umi = Matrix::colSums(qc_mat)
qc_genes = Matrix::colSums(qc_mat > 0)
rm(qc_mat)
ed_genes = gene_type$V1[gene_type$V2 %in% c("lincRNA", "protein_coding", "pseudogene")]   # as in emptydrops.R
ed_umi = Matrix::colSums(rnaseq_data_raw[row.names(rnaseq_data_raw) %in% ed_genes, ])
qc_cells = list(barhop = names(ed_umi)[ed_umi >= as.numeric(barhop_rna)],
                emptydrops = combined_barcodes$V1)
qc_cells$mito = setdiff(qc_cells$emptydrops, high_mt_cells)
mt_mad_thr = tryCatch(attr(stats$high_subsets_Mt_percent, "thresholds")[["higher"]], error = function(e) NA)
rm(ed_umi)

### remove cells with low or too high ATAC fragment counts
counts <- Read10X_h5(paste0("data/sc/cellranger_arc/", sample_id, "/outs/raw_feature_bc_matrix.h5"))
colnames(counts$Peaks) <- sapply(strsplit(colnames(counts$Peaks), split='-', fixed=TRUE), `[`, 1)
atac_count = counts$Peaks
atac_peak_counts = Matrix::colSums(atac_count)   # per barcode, for the QC summary

min_frags_barcodes = colnames(atac_count[,colSums(atac_count) > min_frags])
atac_count = atac_count[,colnames(atac_count) %in% min_frags_barcodes]

# remove barcodes with a number of fragments exceeding 5 MADs
# the median / MAD are computed on the nuclei only (EmptyDrops calls passing the mito filter and the minimum
# above). Computing them on every barcode with > MIN_ATAC_FRAGS counts, as in the paper, includes ~100,000
# near-empty droplets: for PMW941_1 that gave a cut-off of ~4,400 counts and removed a third of the nuclei,
# the best-covered ones, instead of rare outliers.
atac_cs = colSums(atac_count)[colnames(atac_count) %in% colnames(rnaseq_data)]
atac_outlier = isOutlier(atac_cs, nmads = 5, type = "higher")
atac_hi_thr = tryCatch(attr(atac_outlier, "thresholds")[["higher"]], error = function(e) NA)
high_frags_barcodes = names(atac_cs)[which(atac_outlier)]

# plot: ATAC cut sites per nucleus before / after the two ATAC filters (nuclei = EmptyDrops calls passing mito)
atac_plot_df = data.frame(barcode = qc_cells$mito,
                          atac = as.numeric(atac_peak_counts[qc_cells$mito]),
                          umi = as.numeric(qc_umi[qc_cells$mito]))
atac_plot_df$atac[is.na(atac_plot_df$atac)] = 0
atac_plot_df$status = ifelse(atac_plot_df$atac <= min_frags, "removed: too few ATAC cut sites",
                      ifelse(atac_plot_df$barcode %in% high_frags_barcodes, "removed: high-count outlier", "kept"))
atac_plot_df$status = factor(atac_plot_df$status, levels = c("kept", "removed: too few ATAC cut sites", "removed: high-count outlier"))
n_low = sum(atac_plot_df$status == "removed: too few ATAC cut sites")
n_high = sum(atac_plot_df$status == "removed: high-count outlier")
n_kept = sum(atac_plot_df$status == "kept")
atac_cols = c("kept" = "grey40", "removed: too few ATAC cut sites" = "#2c7bb6", "removed: high-count outlier" = "#d7191c")
atac_lines = data.frame(x = c(min_frags, atac_hi_thr),
                        what = c(paste0("MIN_ATAC_FRAGS = ", min_frags),
                                 paste0("high-count cut-off = ", round(atac_hi_thr), " (median + 5 MAD of nuclei)")))
atac_lines = atac_lines[!is.na(atac_lines$x), ]
atac_sub = paste0(nrow(atac_plot_df), " nuclei (EmptyDrops calls after the mito filter): ", n_kept, " kept, ",
                  n_low, " removed with <= ", min_frags, " cut sites, ", n_high, " removed as high-count outliers")
p_atac_hist = ggplot(atac_plot_df, aes(x = atac + 1, fill = status)) +
  geom_histogram(bins = 100, position = "stack") +
  geom_vline(data = atac_lines, aes(xintercept = x + 1, linetype = what), colour = "black") +
  scale_x_log10() + scale_fill_manual(values = atac_cols, drop = FALSE) +
  scale_linetype_manual(values = c("dashed", "dotted")) +
  labs(title = paste0(sample_id, ": ATAC cut sites in Cell Ranger peaks per nucleus"),
       subtitle = atac_sub, x = "ATAC cut sites in peaks + 1 (log scale)", y = "nuclei",
       fill = NULL, linetype = NULL) +
  theme_bw() + theme(legend.position = "bottom", legend.box = "vertical")
p_atac_scatter = ggplot(atac_plot_df[order(atac_plot_df$status), ], aes(x = atac + 1, y = umi + 1, colour = status)) +
  geom_point(size = 0.3, alpha = 0.5) +
  geom_vline(data = atac_lines, aes(xintercept = x + 1, linetype = what), colour = "black") +
  scale_x_log10() + scale_y_log10() + scale_colour_manual(values = atac_cols, drop = FALSE) +
  scale_linetype_manual(values = c("dashed", "dotted")) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = paste0(sample_id, ": ATAC cut sites in peaks vs RNA UMIs per nucleus"),
       subtitle = "RNA UMIs on the genes kept for analysis (lincRNA, protein-coding, pseudogene; no MT genes)",
       x = "ATAC cut sites in peaks + 1 (log scale)", y = "RNA UMIs + 1 (log scale)", colour = NULL, linetype = NULL) +
  theme_bw() + theme(legend.position = "bottom", legend.box = "vertical")
pdf(file = paste0(plot_dir_w_date, "/", sample_id, ".", annotation, ".ATAC_filters.pdf"), width = 10, height = 8)
print(p_atac_hist)
print(p_atac_scatter)
dev.off()
rm(atac_plot_df)
counts = NULL
atac_count = NULL

rnaseq_data <- rnaseq_data[,colnames(rnaseq_data) %in% min_frags_barcodes & colnames(rnaseq_data) %nin% high_frags_barcodes]
message(sample_id, ": ", ncol(rnaseq_data), " barcodes after mito / ATAC filters")
qc_cells$atac_min = intersect(qc_cells$mito, min_frags_barcodes)
qc_cells$atac_high = colnames(rnaseq_data)

combined_data_seurat <- CreateSeuratObject(counts = rnaseq_data, project = sample_id)
combined_data_seurat <- SCTransform(combined_data_seurat, vst.flavor = "v2", method = "glmGamPoi")
combined_data_seurat <- RunPCA(combined_data_seurat, seed.use=42, features = VariableFeatures(object = combined_data_seurat))
combined_data_seurat <- FindNeighbors(combined_data_seurat, dims = 1:50)
combined_data_seurat <- FindClusters(combined_data_seurat, resolution = 2, random.seed = 0)
combined_data_seurat <- RunUMAP(combined_data_seurat, dims = 1:50, seed.use = 42)

pdf(file=paste0(plot_dir_w_date, "/", sample_id, ".", annotation, ".umap_pre_doublets.50PCA.pdf"), width=8, height = 8)
print(DimPlot(combined_data_seurat, reduction = "umap", label=T) +
      labs(title = paste0(sample_id, ": RNA-only UMAP before doublet removal"),
           subtitle = "SCT normalisation, 50 PCA dims"))
dev.off()

# identify doublets
sce <- as.SingleCellExperiment(combined_data_seurat)
set.seed(42)
sce.mam.dbl <- scDblFinder(sce, clusters=colData(sce)$ident)
pdf(file=paste0(plot_dir_w_date, "/", sample_id, ".", annotation, ".seurat.scDblFinder.scores.pdf"), width=8, height = 8)
print(plotUMAP(sce.mam.dbl, colour_by="scDblFinder.score") +
      labs(title = paste0(sample_id, ": scDblFinder doublet score"),
           subtitle = "High score = likely doublet; nuclei classified as doublets are removed"))
dev.off()
print(table(sce.mam.dbl$scDblFinder.class))

singlet_cells = colnames(sce.mam.dbl[,sce.mam.dbl$scDblFinder.class == "singlet"])
rnaseq_data_nodub = rnaseq_data[,colnames(rnaseq_data) %in% singlet_cells]

### additionally, remove cells with extremely high RNA counts
rnaseq_data_nodub = rnaseq_data_nodub[,colSums(rnaseq_data_nodub) < max_rna_umi]
qc_cells$doublets = intersect(qc_cells$atac_high, singlet_cells)
qc_cells$umi_cap = colnames(rnaseq_data_nodub)
n_dbl = sum(sce.mam.dbl$scDblFinder.class == "doublet")
dbl_thr = tryCatch(S4Vectors::metadata(sce.mam.dbl)$scDblFinder.threshold, error = function(e) NULL)
if (is.null(dbl_thr)) dbl_thr = if (n_dbl > 0) min(sce.mam.dbl$scDblFinder.score[sce.mam.dbl$scDblFinder.class == "doublet"]) else NA

combined_data_seurat_nodoub <- CreateSeuratObject(counts = rnaseq_data_nodub, project = sample_id)
combined_data_seurat_nodoub <- SCTransform(combined_data_seurat_nodoub, vst.flavor = "v2", method = "glmGamPoi")
combined_data_seurat_nodoub <- RunPCA(combined_data_seurat_nodoub, seed.use=42)
combined_data_seurat_nodoub <- FindNeighbors(combined_data_seurat_nodoub, dims = 1:50)
combined_data_seurat_nodoub <- FindClusters(combined_data_seurat_nodoub, resolution = 2, random.seed = 0)
combined_data_seurat_nodoub <- RunUMAP(combined_data_seurat_nodoub, dims = 1:50, seed.use = 42)

pdf(file=paste0(plot_dir_w_date, "/", sample_id, ".", annotation, ".umap_pre_soupx.50PCA.pdf"), width=8, height = 8)
print(DimPlot(combined_data_seurat_nodoub, reduction = "umap", label=T) +
      labs(title = paste0(sample_id, ": RNA-only UMAP after doublet removal, before SoupX"),
           subtitle = "SCT normalisation, 50 PCA dims"))
dev.off()

### germline clusters, excluded from the contamination estimate --------------------------------------
# (the paper excluded P-cell clusters found with embryonic P-granule transcripts; for L3 we use
#  the germline panel of the L3 tissue marker consensus)
rnaseq_data_raw_filtered <- rnaseq_data_raw[row.names(rnaseq_data_raw) %in% genes_to_retain,]
germline_genes = intersect(germline_genes_all, rownames(combined_data_seurat_nodoub))
message(sample_id, ": ", length(germline_genes), " of ", length(germline_genes_all), " L3 germline markers found in the matrix")

combined_data_seurat_nodoub = AddModuleScore(combined_data_seurat_nodoub, features = list(germline_genes), nbin = 10, name = "germline_score")

DefaultAssay(combined_data_seurat_nodoub) = "RNA"
germline_avgExp = as.data.frame(AverageExpression(combined_data_seurat_nodoub, slot = "data", assays = "RNA", features = germline_genes, group.by = "seurat_clusters"))
germline_avgExp_log_avg = apply(log10(germline_avgExp + 1), 2, mean)
names(germline_avgExp_log_avg) = levels(combined_data_seurat_nodoub$seurat_clusters)
DefaultAssay(combined_data_seurat_nodoub) = "SCT"

pdf(file = paste0(plot_dir_w_date, "/", sample_id, ".", annotation, ".germline_binarisation.pdf"), width = 8, height = 8, useDingbats = F)
low_high_exp_loc = tryCatch({
  loc <- locmodes(as.numeric(germline_avgExp_log_avg), mod0 = 2, display = F)
  second_derivative_binarization(as.numeric(germline_avgExp_log_avg), loc$cbw$bw, loc$locations, "L3 germline markers (mean per cluster)")
}, error = function(e) {
  warning("germline binarisation failed (", conditionMessage(e), "); no germline cluster excluded from SoupX estimation")
  c(Inf, Inf)
})
dev.off()
# clusters in the high-expression mode (above the second threshold) are excluded
clusters_to_remove = names(germline_avgExp_log_avg)[germline_avgExp_log_avg >= low_high_exp_loc[2]]
clusters_to_remove = clusters_to_remove[!is.na(clusters_to_remove)]
message("germline clusters excluded from contamination estimate: ", paste(clusters_to_remove, collapse = ","))

pdf(file=paste0(plot_dir_w_date, "/", sample_id, ".", annotation, ".germline_clusters.noDub.pdf"), width=16, height = 8)
if (length(clusters_to_remove) > 0) {
  p1 = (DimPlot(combined_data_seurat_nodoub, reduction = "umap", label=T,
                cells.highlight = colnames(combined_data_seurat_nodoub)[combined_data_seurat_nodoub$seurat_clusters %in% clusters_to_remove],
                pt.size = 0.3, sizes.highlight = 0.3) +
        labs(title = paste0(sample_id, ": germline clusters (highlighted)"),
             subtitle = paste0("Clusters ", paste(clusters_to_remove, collapse = ", "),
                               " excluded from the SoupX contamination estimate (kept in the data)"))) &
       NoLegend()
} else {
  p1 = (DimPlot(combined_data_seurat_nodoub, reduction = "umap", label=T, pt.size = 0.3) +
        labs(title = paste0(sample_id, ": no germline cluster identified"),
             subtitle = "Binarisation found no high-expression mode: all clusters used for the SoupX estimate")) &
       NoLegend()
}
p2 = (FeaturePlot(combined_data_seurat_nodoub, reduction = "umap", features = "germline_score1", order = TRUE) +
      labs(title = "L3 germline marker module score",
           subtitle = paste0(length(germline_genes), " germline genes from ", basename(l3_marker_file)))) &
     scale_colour_gradientn(colours = feature_colours)
print(p1 | p2)
dev.off()

### remove barcodes from flagged clusters, then run SoupX
flagged_barcodes = colnames(combined_data_seurat_nodoub)[combined_data_seurat_nodoub$seurat_clusters %in% clusters_to_remove]
rnaseq_data_nodub_filtered = rnaseq_data_nodub[,colnames(rnaseq_data_nodub) %nin% flagged_barcodes]
if (length(clusters_to_remove) > 0) {
  combined_data_seurat_nodoub_filtered = subset(x = combined_data_seurat_nodoub, idents = clusters_to_remove, invert = TRUE)
} else {
  combined_data_seurat_nodoub_filtered = combined_data_seurat_nodoub
}

sc = SoupChannel(rnaseq_data_raw_filtered, rnaseq_data_nodub_filtered, calcSoupProfile=F)
sc = estimateSoup(sc,soupRange=c(as.numeric(barhop_rna) ,as.numeric(lower_rna))  )

### soup profile: which genes make up the ambient RNA ------------------------------------------------
# (MT genes, rRNAs and small RNAs were removed above; see soup_marker_fractions.pdf for those)
soup_prof = sc$soupProfile[order(-sc$soupProfile$est), ]
soup_prof = data.frame(gene = rownames(soup_prof), display_name = disp(rownames(soup_prof)), soup_rank = seq_len(nrow(soup_prof)),
                       est = soup_prof$est, counts = soup_prof$counts,
                       l3_tissue = unname(tissue_of_gene[rownames(soup_prof)]), stringsAsFactors = FALSE)
write.table(soup_prof, file = paste0(seurat_objects_dir, sample_id, ".", annotation, ".soupProfile.tsv"),
            sep = "\t", quote = F, row.names = F)

top_soup = head(soup_prof, 40)
top_soup$label = factor(make.unique(top_soup$display_name), levels = rev(make.unique(top_soup$display_name)))
top_soup$l3_tissue[is.na(top_soup$l3_tissue)] = "not a single-tissue L3 marker"
pdf(file=paste0(plot_dir_w_date, "/", sample_id, ".", annotation, ".soup_profile_top40.pdf"), width=9, height = 10)
print(ggplot(top_soup, aes(x = label, y = 100 * est, fill = l3_tissue)) + geom_col() + coord_flip() +
      labs(title = paste0(sample_id, ": top 40 genes in the ambient RNA (SoupX soup profile)"),
           subtitle = paste0("% of soup UMIs per gene, estimated from droplets with ", barhop_rna, "-", lower_rna,
                             " UMIs.\nColour = L3 tissue for single-tissue markers; full table in seurat_objects/*.soupProfile.tsv"),
           x = NULL, y = "% of soup UMIs", fill = "L3 tissue") +
      theme_classic() + theme(plot.subtitle = element_text(size = 8)))
dev.off()

# tissue-specific soup genes for the change maps: highest-ranking soup gene of each tissue,
# then the top N of those by soup abundance (keeps the set spread over tissues)
soup_tissue = soup_prof[!is.na(soup_prof$l3_tissue) & soup_prof$counts > 0, ]
soup_change_genes = head(soup_tissue[!duplicated(soup_tissue$l3_tissue), "gene"], n_soup_change)
# genes for a global soup score: top soup genes overall
soup_score_genes = head(soup_prof$gene[soup_prof$counts > 0], n_soup_score)
message(sample_id, ": tissue-specific soup genes for change maps: ",
        paste0(disp(soup_change_genes), " (", tissue_of_gene[soup_change_genes], ")", collapse = ", "))

meta    <- combined_data_seurat_nodoub_filtered@meta.data
umap    <- combined_data_seurat_nodoub_filtered@reductions$umap@cell.embeddings
sc  <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
sc  <- setDR(sc, umap)

pdf(file=paste0(plot_dir_w_date, "/",sample_id,".", annotation, ".soupx_contamination_rate.pdf"), width=8, height = 8)
par(mfrow=c(1,1), mar=c(6, 4, 3, 3), oma=c(0,0,2,0))
sc  <- autoEstCont(sc, soupQuantile = soup_quantile)
mtext(paste0(sample_id, ": SoupX ambient RNA contamination estimate (soup quantile = ", soup_quantile, ")"),
      outer = TRUE, line = 0.5, cex = 0.9)
dev.off()

estimated_contamination = unique(sc$metaData$rho)
message(sample_id, ": estimated contamination ", estimated_contamination)

# new SoupChannel including the germline cells, corrected with the estimated contamination
sc_full = SoupChannel(rnaseq_data_raw_filtered, rnaseq_data_nodub, calcSoupProfile=F)
sc_full = estimateSoup(sc_full,soupRange=c(as.numeric(barhop_rna) ,as.numeric(lower_rna))  )

meta    <- combined_data_seurat_nodoub@meta.data
umap    <- combined_data_seurat_nodoub@reductions$umap@cell.embeddings
sc_full  <- setClusters(sc_full, setNames(meta$seurat_clusters, rownames(meta)))
sc_full  <- setDR(sc_full, umap)

sc_full = setContaminationFraction(sc_full, estimated_contamination)
adj.matrix  <- adjustCounts(sc_full, roundToInt = T)

# step-by-step QC summary (stage 5 rows)
qc_write(qc_stage5_rows(
  cells = qc_cells, umi = qc_umi, genes = qc_genes, atac = atac_peak_counts,
  p = list(barhop = barhop_rna, ambient = lower_rna, mt_max = mt_pct_threshold, mt_mad = mt_mad_thr,
           min_frags = min_frags, atac_hi = atac_hi_thr, dbl_thr = dbl_thr, n_dbl = n_dbl, max_umi = max_rna_umi,
           rho = estimated_contamination, soup_quantile = soup_quantile, germline_clusters = clusters_to_remove),
  soup_umi = Matrix::colSums(adj.matrix), soup_genes = Matrix::colSums(adj.matrix > 0)),
  qc_summary_file(sample_id, annotation))

### decontamination change maps for tissue-specific soup genes (replaces embryonic vet-6) ------------
pdf(file=paste0(plot_dir_w_date, "/",sample_id,".", annotation, ".soupx_decontamination.", estimated_contamination, ".soup_marker_change_maps.pdf"), width=16, height = 8)
for (g in soup_change_genes) {
  tryCatch(print(plotChangeMap(sc_full, adj.matrix, geneSet = g, dataType = "binary") +
    labs(title = paste0(sample_id, ": ", disp(g), " (", pretty_tissue(tissue_of_gene[g]), " marker, soup rank ",
                        soup_prof$soup_rank[soup_prof$gene == g], ") detected before / after SoupX"),
         subtitle = paste0("Contamination fraction ", round(estimated_contamination, 3),
                           ". A nucleus is coloured if it has >= 1 ", disp(g), " UMI. After correction, detection outside ",
                           pretty_tissue(tissue_of_gene[g]), " clusters should largely disappear."))),
    error = function(e) warning("plotChangeMap failed for ", g, ": ", conditionMessage(e)))
}
tryCatch(print(plotChangeMap(sc_full, adj.matrix, geneSet = soup_change_genes, dataType = "soupFrac") +
  labs(title = paste0(sample_id, ": fraction of the counts of these genes removed by SoupX, per nucleus"),
       subtitle = paste0(paste(disp(soup_change_genes), collapse = ", "),
                         ". High values (most counts removed) are expected outside the tissues expressing them."))),
  error = function(e) warning("plotChangeMap soupFrac failed: ", conditionMessage(e)))
dev.off()

# re-run clustering with Seurat
combined_data_seurat_soupx <- CreateSeuratObject(counts = adj.matrix, project = paste0(sample_id,"-soupX") )
combined_data_seurat_soupx <- SCTransform(combined_data_seurat_soupx, vst.flavor = "v2", method = "glmGamPoi")
combined_data_seurat_soupx <- RunPCA(combined_data_seurat_soupx, seed.use = 42, features = VariableFeatures(object = combined_data_seurat_soupx), npcs = 100)
combined_data_seurat_soupx <- FindNeighbors(combined_data_seurat_soupx, dims = 1:50)
combined_data_seurat_soupx <- FindClusters(combined_data_seurat_soupx, resolution = 3, random.seed = 0)
combined_data_seurat_soupx <- RunUMAP(combined_data_seurat_soupx, seed.use = 42, dims = 1:50, return.model=TRUE, reduction.name = "umap_SCT")

### QC metadata: soup markers, soup score, germline score ---------------------------------------------
cells_soupx = colnames(combined_data_seurat_soupx)
combined_data_seurat_soupx$mito_pct <- mt_pct_raw[cells_soupx]
combined_data_seurat_soupx$rrna_pct <- rrna_pct_raw[cells_soupx]
combined_data_seurat_soupx$ribo_prot_pct <- ribo_pct_raw[cells_soupx]
sct_data = LayerData(combined_data_seurat_soupx, assay = "SCT", layer = "data")
mean_score = function(genes) { g = intersect(genes, rownames(sct_data)); if (length(g) == 0) rep(0, ncol(sct_data)) else Matrix::colMeans(sct_data[g, , drop = FALSE]) }
combined_data_seurat_soupx$soup_score     <- mean_score(soup_score_genes)
combined_data_seurat_soupx$germline_score <- mean_score(germline_genes)
rm(sct_data)

write.table(combined_data_seurat_soupx[[]], file = paste0(seurat_objects_dir, sample_id, ".", annotation, ".postSoupX.table.txt"), quote = F, sep = "\t", row.names = T, col.names = T)

pdf(file=paste0(plot_dir_w_date, "/",sample_id,".", annotation, ".umap_post_soupx.", estimated_contamination, "_soupX.50PCA.pdf"), width=8, height = 8)
print(DimPlot(combined_data_seurat_soupx, reduction = "umap_SCT", group.by = "seurat_clusters", label = TRUE, label.size = 4) +
      labs(title = paste0(sample_id, ": RNA-only UMAP after SoupX decontamination"),
           subtitle = paste0("Contamination fraction: ", round(estimated_contamination, 3),
                             "  |  SCT normalisation, 50 PCA dims")))
dev.off()

### L3 tissue markers before / after SoupX (replaces embryonic lineage markers) ------------------------
pdf(file=paste0(plot_dir_w_date, "/",sample_id,".", annotation, ".umap_post_soupx.", estimated_contamination, "_soupX.L3_tissue_markers.pdf"),
    width = 5 * n_markers_plot, height = 10)
for (t in l3_tissues) {
  genes_t = l3_single$gene[l3_single$canonical_tissue == t]
  genes_t = head(genes_t[genes_t %in% rownames(combined_data_seurat_soupx) & genes_t %in% rownames(combined_data_seurat_nodoub)], n_markers_plot)
  if (length(genes_t) == 0) { message("no ", t, " markers in the matrix"); next }
  p_pre  = FeaturePlot(combined_data_seurat_nodoub, reduction = "umap", features = genes_t, order = TRUE, ncol = n_markers_plot)
  p_post = FeaturePlot(combined_data_seurat_soupx, reduction = "umap_SCT", features = genes_t, order = TRUE, ncol = n_markers_plot)
  # panel titles: common names (a single-gene FeaturePlot is a plain ggplot, not a patchwork)
  if (inherits(p_pre, "patchwork")) {
    for (i in seq_along(genes_t)) {
      p_pre[[i]]  = p_pre[[i]]  + ggtitle(paste0(disp(genes_t[i]), " (before SoupX)"))
      p_post[[i]] = p_post[[i]] + ggtitle(paste0(disp(genes_t[i]), " (after SoupX)"))
    }
  } else {
    p_pre  = p_pre  + ggtitle(paste0(disp(genes_t[1]), " (before SoupX)"))
    p_post = p_post + ggtitle(paste0(disp(genes_t[1]), " (after SoupX)"))
  }
  p_pre  = p_pre  & scale_colour_gradientn(colours = feature_colours)
  p_post = p_post & scale_colour_gradientn(colours = feature_colours)
  print((p_pre / p_post) +
        plot_annotation(title = paste0(sample_id, ": top L3 ", pretty_tissue(t), " markers"),
                        subtitle = paste0("Top row: before SoupX (UMAP of doublet-filtered nuclei).  Bottom row: after SoupX (post-SoupX UMAP).\n",
                                          "Markers ranked by the L3 consensus (", basename(l3_marker_file), "). ",
                                          "Expression should become more restricted to one region after SoupX.")))
}
dev.off()

pdf(file=paste0(plot_dir_w_date, "/",sample_id,".", annotation, ".umap_post_soupx.", estimated_contamination, "_soupX.qc_bad_clusters.pdf"), width=20, height = 10)
print((FeaturePlot(combined_data_seurat_soupx, reduction = "umap_SCT", order = T, ncol = 4,
                   features = c("mito_pct", "rrna_pct", "ribo_prot_pct", "soup_score", "nFeature_RNA", "nCount_RNA", "germline_score")) &
       scale_colour_gradientn(colours = feature_colours)) +
      plot_annotation(title = paste0(sample_id, ": QC metrics for identifying low-quality or soup-dominated clusters (after SoupX)"),
                      subtitle = paste0("mito / rRNA / ribo_prot = % of raw UMIs (mito is cytoplasmic, i.e. mostly soup-derived in nuclei); ",
                                        "soup_score = mean expression of the top ", n_soup_score, " soup genes (post-SoupX);\n",
                                        "germline_score = mean expression of L3 germline markers. Clusters with high soup markers and ",
                                        "low nFeature_RNA are likely damaged nuclei or empty droplets.")))
dev.off()

# add ATAC data quantified over the bulk scATAC peak set
merged_peaks = import.bed(merged_peaks_file)
fragpath_raw <- paste0("data/sc/cellranger_arc/", sample_id, "/outs/atac_fragments_no_dash.tsv.gz")
# cellranger-arc produces a 6-column fragment file (chrom/start/end/barcode/count/strand).
# Signac requires exactly 5 columns, so strip the strand column if present.
fragpath <- sub("\\.tsv\\.gz$", "_5col.tsv.gz", fragpath_raw)
if (!file.exists(fragpath)) {
  message(sample_id, ": stripping strand column from fragment file for Signac compatibility...")
  tmp_tsv <- tempfile(pattern = paste0(sample_id, "_frags5col_"), tmpdir = "/tmp", fileext = ".tsv")
  ret <- system(paste0("zcat ", fragpath_raw, " | cut -f1-5 > ", tmp_tsv))
  if (ret != 0) stop("cut -f1-5 failed when creating 5-column fragment file")
  Rsamtools::bgzip(tmp_tsv, dest = fragpath, overwrite = TRUE)
  Rsamtools::indexTabix(fragpath, format = "bed")
  file.remove(tmp_tsv)
  message(sample_id, ": 5-column fragment file written to ", fragpath)
}

fragments_sample <- CreateFragmentObject(fragpath, cells = colnames(combined_data_seurat_soupx))
merged_peaks_matrix_filtered <- FeatureMatrix(fragments_sample, sep = c("-", "-"), features=merged_peaks, cells = colnames(combined_data_seurat_soupx))

combined_data_seurat_soupx[["ATAC"]] <- CreateChromatinAssay(
  counts = merged_peaks_matrix_filtered,
  sep = c(":", "-"),
  fragments = fragments_sample)

DefaultAssay(combined_data_seurat_soupx) <- "SCT"
pdf(file=paste0(plot_dir_w_date, "/",sample_id,".", annotation, ".umap_post_soupx.", estimated_contamination, "_soupX.ATAC_features_counts.pdf"), width=24, height = 8)
p1 = DimPlot(combined_data_seurat_soupx, reduction = "umap_SCT", group.by = "seurat_clusters", label = TRUE, label.size = 4) +
     labs(title = "RNA clusters (SCT)")
p2 = FeaturePlot(combined_data_seurat_soupx, features = "nFeature_ATAC", min.cutoff = 0, max.cutoff = 3000, order=T) & scale_colour_gradientn(colours = feature_colours)
p3 = FeaturePlot(combined_data_seurat_soupx, features = "nCount_ATAC", min.cutoff = 0, max.cutoff = 3000, order=T) & scale_colour_gradientn(colours = feature_colours)
print(p1 + p2 + p3 + plot_annotation(
  title = paste0(sample_id, ": ATAC chromatin accessibility QC – peaks and counts per nucleus (capped at 3 000)")))
dev.off()

ATAC_counts <- CountFragments(fragments = fragpath, cells = colnames(combined_data_seurat_soupx))
DefaultAssay(combined_data_seurat_soupx) <- "ATAC"
combined_data_seurat_soupx$ATAC_count <- ATAC_counts$reads_count[match(colnames(combined_data_seurat_soupx), ATAC_counts$CB)]
combined_data_seurat_soupx <- FRiP(combined_data_seurat_soupx, assay = "ATAC", total.fragments = "ATAC_count", col.name = "FRiP", verbose = TRUE)

pdf(file=paste0(plot_dir_w_date, "/",sample_id,".", annotation, ".umap_post_soupx.", estimated_contamination, "_soupX.ATAC_features_FRiP.pdf"), width=24, height = 8)
p1 = DimPlot(combined_data_seurat_soupx, reduction = "umap_SCT", group.by = "seurat_clusters", label = TRUE, label.size = 4) +
     labs(title = "RNA clusters (SCT)")
p2 = FeaturePlot(combined_data_seurat_soupx, features = "nFeature_ATAC", min.cutoff = 0, max.cutoff = 3000, order=T) & scale_colour_gradientn(colours = feature_colours)
p3 = FeaturePlot(combined_data_seurat_soupx, features = "FRiP", min.cutoff = 0, max.cutoff = 0.3, order=T) & scale_colour_gradientn(colours = feature_colours)
print(p1 + p2 + p3 + plot_annotation(
  title = paste0(sample_id, ": ATAC QC – Fraction of Reads in Peaks (FRiP, capped at 0.30)"),
  subtitle = "FRiP reflects ATAC signal-to-noise; low values (< 0.15) suggest poor-quality ATAC data"))
dev.off()

# ATAC-seq clustering
DefaultAssay(combined_data_seurat_soupx) = "ATAC"
combined_data_seurat_soupx <- RunTFIDF(combined_data_seurat_soupx)
combined_data_seurat_soupx <- FindTopFeatures(combined_data_seurat_soupx, min.cutoff = 'q0')
combined_data_seurat_soupx <- RunSVD(combined_data_seurat_soupx)
combined_data_seurat_soupx <- RunUMAP(object = combined_data_seurat_soupx, reduction = 'lsi', dims = 2:30, reduction.name = "umap_ATAC")
combined_data_seurat_soupx <- FindNeighbors(object = combined_data_seurat_soupx, reduction = 'lsi', dims = 2:30)
combined_data_seurat_soupx <- FindClusters(object = combined_data_seurat_soupx, verbose = T, algorithm = 3, resolution = 3)

pdf(file=paste0(plot_dir_w_date, "/",sample_id,".", annotation, ".umap_ATAC_MACS_peaks.pdf"), width=8, height = 8)
print(DimPlot(object = combined_data_seurat_soupx, label = TRUE) +
      labs(title = paste0(sample_id, ": ATAC-only UMAP (bulk peak set, LSI dims 2-30)"),
           subtitle = "Independent ATAC clustering; compare with the RNA-based clusters"))
dev.off()

saveRDS(combined_data_seurat_soupx, file = paste0(seurat_objects_dir, sample_id, ".", annotation, ".postSoupX.ATAC_MACS.rds"))
message(sample_id, ": ", ncol(combined_data_seurat_soupx), " nuclei saved")
