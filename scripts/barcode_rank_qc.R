# Helper (new): UMI-per-barcode plots to choose the 'barhop' and 'ambient' thresholds of each sample.
# In the paper these were set manually from the UMI distribution (see emptydrops.20231106.R).
#   barhop  = upper end of the barcode-hopping peak (barcodes below are discarded before EmptyDrops)
#   ambient = EmptyDrops 'lower' bound (barcodes below are used to model ambient RNA)
# Paper values: exp024 100/200, exp025 100/400, exp027 50/200, exp031 100/400, exp032 100/400, exp042 20/150, exp043 20/150
suppressPackageStartupMessages({ library(Matrix); library(DropletUtils) })
source("scripts/pipeline_utils.R")

args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
wb_annotation <- args[2]
canon <- canon_prefix()
plot_dir <- make_plot_dir("barcode_rank_qc")

m <- Seurat::Read10X(paste0("data/sc/rnaseq/", sample_id, "/star_", wb_annotation, "/star_", wb_annotation, "_Solo.out/GeneFull/raw/um_reads"))
genes_type <- read.table(paste0(canon, ".gene_type.txt"))
gene_name <- read.table(paste0(canon, ".gene_name.txt"))
keep <- gene_name$V2[gene_name$V1 %in% genes_type$V1[genes_type$V2 %in% c("lincRNA", "protein_coding", "pseudogene")]]
m <- m[rownames(m) %in% keep, ]
counts <- Matrix::colSums(m)
counts <- counts[counts > 0]

br <- barcodeRanks(m, lower = 10)
knee <- metadata(br)$knee
inflection <- metadata(br)$inflection

pdf(file = file.path(plot_dir, paste0(sample_id, ".", wb_annotation, ".umi_distribution.pdf")), width = 10, height = 6)
hist(counts, breaks = seq(0, ceiling(max(counts)) + 1, by = 1), xlim = c(0, 500), ylim = c(0, 10000), las = 1,
     main = paste(sample_id, "- UMIs per barcode (0-500)"), xlab = "UMIs")
abline(v = seq(0, 500, 50), col = "grey80", lty = 3)
hist(log10(counts), breaks = 200, las = 1, main = paste(sample_id, "- log10 UMIs per barcode"), xlab = "log10(UMIs)")
o <- order(br$rank)
plot(br$rank[o], br$total[o], log = "xy", xlab = "barcode rank", ylab = "UMIs", main = paste(sample_id, "- barcode rank plot"), pch = 16, cex = 0.3)
abline(h = knee, col = "dodgerblue", lty = 2); abline(h = inflection, col = "forestgreen", lty = 2)
abline(h = c(20, 50, 100, 150, 200, 400), col = "grey70", lty = 3)
legend("bottomleft", legend = c(paste("knee =", round(knee)), paste("inflection =", round(inflection))),
       col = c("dodgerblue", "forestgreen"), lty = 2, bty = "n")
dev.off()

cat(sprintf("sample\t%s\nbarcodes_with_umis\t%d\nknee\t%.0f\ninflection\t%.0f\nbarcodes_above_inflection\t%d\n",
            sample_id, length(counts), knee, inflection, sum(counts >= inflection)),
    file = file.path(plot_dir, paste0(sample_id, ".", wb_annotation, ".summary.txt")))
