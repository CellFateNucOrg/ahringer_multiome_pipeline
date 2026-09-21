# Per-sample QC: mito / ATAC fragment filters, scDblFinder, SoupX, ATAC assay on the bulk peak set
# (adapted from seurat_soupx.20240207.R; commented-out exploratory code removed)
# changes: thresholds from samples.tsv / config, plots in plots/seurat_soupx/<sample>,
#          ATAC_count matched by barcode, fallback when P-cell binarisation is not possible
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

ED_barcodes = paste0("data/sc/rnaseq/", sample_id, "/ED_RNA_", annotation, "/ED_RNA_barcodes.txt")
set.seed(42)

plot_dir_w_date = make_plot_dir(file.path("seurat_soupx", sample_id))
seurat_objects_dir = "seurat_objects/"
dir.create(seurat_objects_dir, showWarnings = FALSE)

raw_matrix_fran = paste0("data/sc/rnaseq/",sample_id, "/star_", annotation, "/star_", annotation, "_Solo.out/GeneFull/raw/um_reads")
rnaseq_data_raw <- Seurat::Read10X(raw_matrix_fran)

# read in the barcodes obtained by EmptyDrops and subset the object
combined_barcodes <- read.table(ED_barcodes)
rnaseq_data <- rnaseq_data_raw[,colnames(rnaseq_data_raw) %in% combined_barcodes$V1]
combined_data_seurat <- CreateSeuratObject(counts = rnaseq_data, project = sample_id )

### inspect mitochondrial contamination
mt_genes <- read.table(paste0(canon, ".filtered.MtDNA_genes.genes"))
gene_type = read.table(paste0(canon, ".gene_name_gene_type.txt"))

mt_genes_expression_pct = colSums(rnaseq_data_raw[row.names(rnaseq_data_raw) %in% mt_genes$V1,])/colSums(rnaseq_data_raw)

sce.rnaseq_data <- as.SingleCellExperiment(combined_data_seurat)
is.mito <- rownames(sce.rnaseq_data) %in% mt_genes$V1
sce.rnaseq_data <- addPerCellQC(sce.rnaseq_data, subsets=list(Mt=is.mito))
print(summary(sce.rnaseq_data$subsets_Mt_percent == 0))

stats <- quickPerCellQC(colData(sce.rnaseq_data), sub.fields="subsets_Mt_percent")

pdf(file=paste0(plot_dir_w_date, "/", sample_id, ".", annotation, ".mitochondrial_pct_all_barcodes.pdf"), width=8, height = 8)
print(plotColData(sce.rnaseq_data, y="subsets_Mt_percent", colour_by=I(stats$high_subsets_Mt_percent)))
dev.off()

high_mt_cells = colnames(sce.rnaseq_data[,sce.rnaseq_data$subsets_Mt_percent >= mt_pct_threshold])

# remove high-Mito cells
rnaseq_data <- rnaseq_data[,colnames(rnaseq_data) %nin% high_mt_cells]

# remove MT genes, rRNAs and small RNAs from matrix
genes_to_retain = gene_type$V1[gene_type$V2 %in% c("lincRNA", "protein_coding", "pseudogene") & gene_type$V1 %nin% mt_genes$V1]
rnaseq_data <- rnaseq_data[row.names(rnaseq_data) %in% genes_to_retain,]

### remove cells with low or too high ATAC fragment counts
counts <- Read10X_h5(paste0("data/sc/cellranger_arc/", sample_id, "/outs/raw_feature_bc_matrix.h5"))
colnames(counts$Peaks) <- sapply(strsplit(colnames(counts$Peaks), split='-', fixed=TRUE), `[`, 1)
atac_count = counts$Peaks

min_frags_barcodes = colnames(atac_count[,colSums(atac_count) > min_frags])
atac_count = atac_count[,colnames(atac_count) %in% min_frags_barcodes]

# remove barcodes with a number of fragments exceeding 5 MADs
high_frags_barcodes = names(colSums(atac_count)[which(isOutlier(colSums(atac_count), nmads = 5))])
counts = NULL
atac_count = NULL

rnaseq_data <- rnaseq_data[,colnames(rnaseq_data) %in% min_frags_barcodes & colnames(rnaseq_data) %nin% high_frags_barcodes]
message(sample_id, ": ", ncol(rnaseq_data), " barcodes after mito / ATAC filters")

combined_data_seurat <- CreateSeuratObject(counts = rnaseq_data, project = sample_id)
combined_data_seurat <- SCTransform(combined_data_seurat, vst.flavor = "v2", method = "glmGamPoi")
combined_data_seurat <- RunPCA(combined_data_seurat, seed.use=42, features = VariableFeatures(object = combined_data_seurat))
combined_data_seurat <- FindNeighbors(combined_data_seurat, dims = 1:50)
combined_data_seurat <- FindClusters(combined_data_seurat, resolution = 2, random.seed = 0)
combined_data_seurat <- RunUMAP(combined_data_seurat, dims = 1:50, seed.use = 42)

pdf(file=paste0(plot_dir_w_date, "/", sample_id, ".", annotation, ".umap_pre_doublets.50PCA.pdf"), width=8, height = 8)
print(DimPlot(combined_data_seurat, reduction = "umap", label=T))
dev.off()

# identify doublets
sce <- as.SingleCellExperiment(combined_data_seurat)
set.seed(42)
sce.mam.dbl <- scDblFinder(sce, clusters=colData(sce)$ident)
pdf(file=paste0(plot_dir_w_date, "/", sample_id, ".", annotation, ".seurat.scDblFinder.scores.pdf"), width=8, height = 8)
print(plotUMAP(sce.mam.dbl, colour_by="scDblFinder.score"))
dev.off()
print(table(sce.mam.dbl$scDblFinder.class))

singlet_cells = colnames(sce.mam.dbl[,sce.mam.dbl$scDblFinder.class == "singlet"])
rnaseq_data_nodub = rnaseq_data[,colnames(rnaseq_data) %in% singlet_cells]

### additionally, remove cells with extremely high RNA counts
rnaseq_data_nodub = rnaseq_data_nodub[,colSums(rnaseq_data_nodub) < max_rna_umi]

combined_data_seurat_nodoub <- CreateSeuratObject(counts = rnaseq_data_nodub, project = sample_id)
combined_data_seurat_nodoub <- SCTransform(combined_data_seurat_nodoub, vst.flavor = "v2", method = "glmGamPoi")
combined_data_seurat_nodoub <- RunPCA(combined_data_seurat_nodoub, seed.use=42)
combined_data_seurat_nodoub <- FindNeighbors(combined_data_seurat_nodoub, dims = 1:50)
combined_data_seurat_nodoub <- FindClusters(combined_data_seurat_nodoub, resolution = 2, random.seed = 0)
combined_data_seurat_nodoub <- RunUMAP(combined_data_seurat_nodoub, dims = 1:50, seed.use = 42)

pdf(file=paste0(plot_dir_w_date, "/", sample_id, ".", annotation, ".umap_pre_soupx.50PCA.pdf"), width=8, height = 8)
print(DimPlot(combined_data_seurat_nodoub, reduction = "umap", label=T))
dev.off()

# remove P-cell clusters before estimating contamination
# get average expression per cluster of P-granule genes
rnaseq_data_raw_filtered <- rnaseq_data_raw[row.names(rnaseq_data_raw) %in% genes_to_retain,]
P_cells = read.table("data/external_data/P_granule_transcripts.sorted.txt")

combined_data_seurat_nodoub = AddModuleScore(combined_data_seurat_nodoub, features = list(P_cells$V1), nbin = 10)

DefaultAssay(combined_data_seurat_nodoub) = "RNA"
combined_data_seurat_nodoub_avgExp = as.data.frame(AverageExpression(combined_data_seurat_nodoub, slot = "data", assays = "RNA", features = P_cells$V1, group.by = "seurat_clusters"))
combined_data_seurat_nodoub_avgExp_log = log10(combined_data_seurat_nodoub_avgExp + 1)
combined_data_seurat_nodoub_avgExp_log_avg = apply(combined_data_seurat_nodoub_avgExp_log, 2, mean)
combined_data_seurat_nodoub_avgExp_quant = combined_data_seurat_nodoub_avgExp_log_avg

pdf(file = paste0(plot_dir_w_date, "/", sample_id, ".", annotation, ".P_cell_binarisation.pdf"), width = 8, height = 8, useDingbats = F)
low_high_exp_loc = tryCatch({
  loc <- locmodes(as.numeric(combined_data_seurat_nodoub_avgExp_log_avg), mod0 = 2, display = F)
  second_derivative_binarization(as.numeric(combined_data_seurat_nodoub_avgExp_log_avg), loc$cbw$bw, loc$locations, "P_granule_genes")
}, error = function(e) {
  warning("P-granule binarisation failed (", conditionMessage(e), "); no P-cell cluster excluded from SoupX estimation")
  c(Inf, Inf)
})
dev.off()
combined_data_seurat_nodoub_avgExp_quant[which(as.numeric(combined_data_seurat_nodoub_avgExp_log_avg) < low_high_exp_loc[1])] <- 0
combined_data_seurat_nodoub_avgExp_quant[which(as.numeric(combined_data_seurat_nodoub_avgExp_log_avg) >= low_high_exp_loc[1] & as.numeric(combined_data_seurat_nodoub_avgExp_log_avg) < low_high_exp_loc[2])] <- 1
combined_data_seurat_nodoub_avgExp_quant[which(as.numeric(combined_data_seurat_nodoub_avgExp_log_avg) >= low_high_exp_loc[2])] <- 2

cl_n = max(as.numeric(combined_data_seurat_nodoub$seurat_clusters))-1
names(combined_data_seurat_nodoub_avgExp_quant) = c(0:cl_n)

clusters_to_remove = names(combined_data_seurat_nodoub_avgExp_quant[combined_data_seurat_nodoub_avgExp_quant == 2])
message("P-cell clusters excluded from contamination estimate: ", paste(clusters_to_remove, collapse = ","))

pdf(file=paste0(plot_dir_w_date, "/", sample_id, ".", annotation, ".seurat_P_cell_clusters.noDub.pdf"), width=16, height = 8)
p1 = DimPlot(combined_data_seurat_nodoub, reduction = "umap", label=T, cells.highlight = colnames(combined_data_seurat_nodoub[,combined_data_seurat_nodoub$seurat_clusters %in% clusters_to_remove]), pt.size = 0.3, sizes.highlight = 0.3) & NoLegend()
p2 = FeaturePlot(combined_data_seurat_nodoub, reduction = "umap", features = "Cluster1")  & scale_colour_gradientn(colours = feature_colours)
print(p1 | p2)
dev.off()

### remove barcodes from flagged clusters, then run SoupX
flagged_barcodes = colnames(combined_data_seurat_nodoub[,combined_data_seurat_nodoub$seurat_clusters %in% clusters_to_remove])
rnaseq_data_nodub_filtered = rnaseq_data_nodub[,colnames(rnaseq_data_nodub) %nin% flagged_barcodes]
if (length(clusters_to_remove) > 0) {
  combined_data_seurat_nodoub_filtered = subset(x = combined_data_seurat_nodoub, idents = clusters_to_remove, invert = TRUE)
} else {
  combined_data_seurat_nodoub_filtered = combined_data_seurat_nodoub
}

sc = SoupChannel(rnaseq_data_raw_filtered, rnaseq_data_nodub_filtered, calcSoupProfile=F)
sc = estimateSoup(sc,soupRange=c(as.numeric(barhop_rna) ,as.numeric(lower_rna))  )

meta    <- combined_data_seurat_nodoub_filtered@meta.data
umap    <- combined_data_seurat_nodoub_filtered@reductions$umap@cell.embeddings
sc  <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
sc  <- setDR(sc, umap)

pdf(file=paste0(plot_dir_w_date, "/",sample_id,".", annotation, ".soupx_contamination_rate.pdf"), width=8, height = 8)
par(mfrow=c(1,1), mar=c(6, 4, 3, 3), oma=c(0,0,0,0))
sc  <- autoEstCont(sc, soupQuantile = soup_quantile)
dev.off()

estimated_contamination = unique(sc$metaData$rho)
message(sample_id, ": estimated contamination ", estimated_contamination)

# new SoupChannel including the P cells, corrected with the estimated contamination
sc_full = SoupChannel(rnaseq_data_raw_filtered, rnaseq_data_nodub, calcSoupProfile=F)
sc_full = estimateSoup(sc_full,soupRange=c(as.numeric(barhop_rna) ,as.numeric(lower_rna))  )

meta    <- combined_data_seurat_nodoub@meta.data
umap    <- combined_data_seurat_nodoub@reductions$umap@cell.embeddings
sc_full  <- setClusters(sc_full, setNames(meta$seurat_clusters, rownames(meta)))
sc_full  <- setDR(sc_full, umap)

sc_full = setContaminationFraction(sc_full, estimated_contamination)
adj.matrix  <- adjustCounts(sc_full, roundToInt = T)

pdf(file=paste0(plot_dir_w_date, "/",sample_id,".", annotation, ".soupx_decontamination.", estimated_contamination, ".vet-6.binary.pdf"), width=16, height = 8)
print(plotChangeMap(sc, adj.matrix, geneSet = c("vet-6"), dataType = "binary"))
dev.off()

# re-run clustering with Seurat
combined_data_seurat_soupx <- CreateSeuratObject(counts = adj.matrix, project = paste0(sample_id,"-soupX") )
combined_data_seurat_soupx <- SCTransform(combined_data_seurat_soupx, vst.flavor = "v2", method = "glmGamPoi")
combined_data_seurat_soupx <- RunPCA(combined_data_seurat_soupx, seed.use = 42, features = VariableFeatures(object = combined_data_seurat_soupx), npcs = 100)
combined_data_seurat_soupx <- FindNeighbors(combined_data_seurat_soupx, dims = 1:50)
combined_data_seurat_soupx <- FindClusters(combined_data_seurat_soupx, resolution = 3, random.seed = 0)
combined_data_seurat_soupx <- RunUMAP(combined_data_seurat_soupx, seed.use = 42, dims = 1:50, return.model=TRUE, reduction.name = "umap_SCT")

write.table(combined_data_seurat_soupx[[]], file = paste0(seurat_objects_dir, sample_id, ".", annotation, ".postSoupX.table.txt"), quote = F, sep = "\t", row.names = T, col.names = T)

pdf(file=paste0(plot_dir_w_date, "/",sample_id,".", annotation, ".umap_post_soupx.", estimated_contamination, "_soupX.50PCA.pdf"), width=8, height = 8)
print(DimPlot(combined_data_seurat_soupx, reduction = "umap_SCT", group.by = "seurat_clusters", label = TRUE, label.size = 4))
dev.off()

pdf(file=paste0(plot_dir_w_date, "/",sample_id,".", annotation, ".umap_post_soupx.", estimated_contamination, "_soupX.markers.pdf"), width=16, height = 16)
print(FeaturePlot(combined_data_seurat_soupx, order = T, features = c("med-1", "end-3", "tbx-35", "elt-7", "pes-10", "lsl-1")) & scale_colour_gradientn(colours = feature_colours))
print(FeaturePlot(combined_data_seurat_soupx, order = T, features = c("end-3", "end-1", "elt-7", "elt-2", "ceh-51", "lsl-1")) & scale_colour_gradientn(colours = feature_colours))
print(FeaturePlot(combined_data_seurat_soupx, order = T, features = c("tbx-33", "tbx-38", "tbx-40", "sdz-20", "pes-10", "vet-2")) & scale_colour_gradientn(colours = feature_colours))
dev.off()

combined_data_seurat_soupx$mito_pct <-  mt_genes_expression_pct[colnames(combined_data_seurat_soupx)]

pdf(file=paste0(plot_dir_w_date, "/",sample_id,".", annotation, ".umap_post_soupx.", estimated_contamination, "_soupX.qc_bad_clusters.pdf"), width=16, height = 16)
print(FeaturePlot(combined_data_seurat_soupx, order = T, features = c("mito_pct", "nFeature_RNA")) & scale_colour_gradientn(colours = feature_colours))
dev.off()

# add ATAC data quantified over the bulk scATAC peak set
merged_peaks = import.bed(merged_peaks_file)
fragpath <- paste0("data/sc/cellranger_arc/", sample_id, "/outs/atac_fragments_no_dash.tsv.gz")

fragments_sample <- CreateFragmentObject(fragpath, cells = colnames(combined_data_seurat_soupx))
merged_peaks_matrix_filtered <- FeatureMatrix(fragments_sample, sep = c("-", "-"), features=merged_peaks, cells = colnames(combined_data_seurat_soupx))

combined_data_seurat_soupx[["ATAC"]] <- CreateChromatinAssay(
  counts = merged_peaks_matrix_filtered,
  sep = c(":", "-"),
  fragments = fragments_sample)

DefaultAssay(combined_data_seurat_soupx) <- "SCT"
pdf(file=paste0(plot_dir_w_date, "/",sample_id,".", annotation, ".umap_post_soupx.", estimated_contamination, "_soupX.ATAC_features_counts.pdf"), width=24, height = 8)
p1 = DimPlot(combined_data_seurat_soupx, reduction = "umap_SCT", group.by = "seurat_clusters", label = TRUE, label.size = 4)
p2 = FeaturePlot(combined_data_seurat_soupx, features = "nFeature_ATAC", min.cutoff = 0, max.cutoff = 3000, order=T) & scale_colour_gradientn(colours = feature_colours)
p3 = FeaturePlot(combined_data_seurat_soupx, features = "nCount_ATAC", min.cutoff = 0, max.cutoff = 3000, order=T) & scale_colour_gradientn(colours = feature_colours)
print(p1 + p2 + p3)
dev.off()

ATAC_counts <- CountFragments(fragments = fragpath, cells = colnames(combined_data_seurat_soupx))
DefaultAssay(combined_data_seurat_soupx) <- "ATAC"
combined_data_seurat_soupx$ATAC_count <- ATAC_counts$reads_count[match(colnames(combined_data_seurat_soupx), ATAC_counts$CB)]
combined_data_seurat_soupx <- FRiP(combined_data_seurat_soupx, assay = "ATAC", total.fragments = "ATAC_count", col.name = "FRiP", verbose = TRUE)

pdf(file=paste0(plot_dir_w_date, "/",sample_id,".", annotation, ".umap_post_soupx.", estimated_contamination, "_soupX.ATAC_features_FRiP.pdf"), width=24, height = 8)
p1 = DimPlot(combined_data_seurat_soupx, reduction = "umap_SCT", group.by = "seurat_clusters", label = TRUE, label.size = 4)
p2 = FeaturePlot(combined_data_seurat_soupx, features = "nFeature_ATAC", min.cutoff = 0, max.cutoff = 3000, order=T) & scale_colour_gradientn(colours = feature_colours)
p3 = FeaturePlot(combined_data_seurat_soupx, features = "FRiP", min.cutoff = 0, max.cutoff = 0.3, order=T) & scale_colour_gradientn(colours = feature_colours)
print(p1 + p2 + p3)
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
print(DimPlot(object = combined_data_seurat_soupx, label = TRUE))
dev.off()

saveRDS(combined_data_seurat_soupx, file = paste0(seurat_objects_dir, sample_id, ".", annotation, ".postSoupX.ATAC_MACS.rds"))
message(sample_id, ": ", ncol(combined_data_seurat_soupx), " nuclei saved")
