# Helpers shared by the adapted R scripts (not part of the original paper code)
suppressPackageStartupMessages(library(Hmisc))  # %nin%

read_samples <- function(path = Sys.getenv("SAMPLES_TSV")) {
  if (path == "") stop("SAMPLES_TSV environment variable is not set")
  s <- read.delim(path, comment.char = "#", colClasses = "character", strip.white = TRUE, check.names = FALSE)
  s[!is.na(s$sample_id) & s$sample_id != "", , drop = FALSE]
}

sample_value <- function(samples, id, col) {
  v <- samples[samples$sample_id == id, col]
  if (length(v) != 1) stop(sprintf("sample '%s' / column '%s' not found in samples.tsv", id, col))
  v
}

env_num <- function(name, default) { v <- Sys.getenv(name, ""); if (v == "") default else as.numeric(v) }
env_chr <- function(name, default = "") { v <- Sys.getenv(name, ""); if (v == "") default else v }

canon_prefix <- function() {
  paste0("species/elegans/gene_annotation/c_elegans.PRJNA13758.", env_chr("WS", "WS285"), ".canonical_geneset")
}

make_plot_dir <- function(name) {
  d <- file.path("plots", name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# merge() that also works with a single object
merge_seurat_list <- function(objs, add.cell.ids = NULL, ...) {
  if (length(objs) == 1) {
    o <- objs[[1]]
    if (!is.null(add.cell.ids)) o <- RenameCells(o, add.cell.id = add.cell.ids[1])
    return(o)
  }
  merge(objs[[1]], y = objs[-1], add.cell.ids = add.cell.ids, ...)
}

# total fragments per cell, matched to the object's cell order
# (the paper assigned CountFragments() output by position, which is not guaranteed to follow cell order)
total_fragments <- function(obj, assay) {
  paths <- lapply(Signac::Fragments(obj[[assay]]), function(f) Signac::GetFragmentData(f, slot = "path"))
  tf <- Signac::CountFragments(paths, cells = colnames(obj), verbose = FALSE)
  tf$frequency_count[match(colnames(obj), tf$CB)]
}

# L3 key tissue markers - data-driven selection across 6 studies
# Ranked by: (1) number of studies with strict specificity (other_max==0, tissue_max>=1)
#            (2) total study support, (3) preferring named genes
# Studies used: shared_min2studies, cengen_L1/L4/adult, cao2017_L2, ghaddar2023
key_lineage_markers <- c(
  "inos-1",  "aman-1",  "clec-145",   # Coelomocytes  (2 strict studies)
  "cpr-3",   "klo-1",   "nspc-1",     # Excretory     (1 strict study; only in shared_min2)
  "pgl-1",   "pid-1",   "sip-1",      # Germline      (4 strict studies: cengen_L1/L4/adult + shared_min2)
  "glt-1",   "grl-12",  "cka-1",      # Glia          (2 strict studies)
  "inx-8",   "inx-9",   "C47D12.5",   # Gonad         (2 strict studies)
  "col-107", "dpy-13",  "ram-2",      # Hypodermis    (2 strict studies)
  "asp-1",   "act-5",   "F57F5.1",    # Intestine     (5 strict studies: all datasets)
  "ace-1",   "acs-16",  "bgal-1",     # Muscle        (1 strict study; shared_min2)
  "acc-1",   "acr-2",   "acr-17",     # Neurons       (2 strict studies)
  "marg-1",  "abf-6",   "abu-14"      # Pharynx       (2 strict studies)
)
feature_colours <- c("grey90", "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred")

# plot_key_markers: Assay5-safe version using LayerData + ggplot2 (no FeaturePlot)
plot_key_markers <- function(obj, reduction, file) {
  library(ggplot2)
  library(SeuratObject)
  emb <- as.data.frame(Embeddings(obj, reduction = reduction))
  colnames(emb) <- c("DIM_1", "DIM_2")
  expr <- LayerData(obj, assay = "RNA", layer = "data")
  genes_use <- intersect(key_lineage_markers, rownames(expr))
  if (length(genes_use) == 0) {
    message("plot_key_markers: no key markers found in RNA assay, skipping")
    return(invisible(NULL))
  }
  pdf(file = file, width = 8, height = 8)
  for (g in genes_use) {
    df <- emb
    df$expr <- as.numeric(expr[g, , drop = FALSE])
    df <- df[order(df$expr), ]
    p <- ggplot(df, aes(x = DIM_1, y = DIM_2, colour = expr)) +
      geom_point(size = 0.3) +
      scale_colour_gradientn(colours = feature_colours, name = g) +
      ggtitle(g) + theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
  dev.off()
}

write_cells_per_group <- function(obj, column, file) {
  groups <- unique(as.character(obj@meta.data[[column]]))
  num <- suppressWarnings(as.numeric(groups))
  groups <- if (all(!is.na(num))) groups[order(num)] else sort(groups)
  tab <- data.frame(row.names = groups, n_cells = as.integer(table(obj@meta.data[[column]])[groups]))
  write.table(tab, file = file, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
}

# TRUE for WormBase common (CGC-style) gene names such as "amn-1", "cpr-9", "bus-8B";
# FALSE for sequence names such as "Y75B7AL.2" or "F57F5.1". Used to prefer common names in plots.
is_common_name <- function(g) grepl("^[A-Za-z]{2,6}-[0-9]", g)

# ---- step-by-step QC summary ------------------------------------------------------------------------
# One TSV per sample, qc_summary/<sample>.<annotation>.qc_steps.tsv, one row per cleaning step:
# number of barcodes / nuclei left, how many the step removed, RNA UMIs, genes and ATAC counts per nucleus
# (median and interquartile range) and the thresholds / parameters the step used.
# Stage 5 (seurat_soupx.R) writes the file; stage 7 (batch_enriched_clusters_removal.R) appends its rows.
QC_BASE_COLS <- c("stage", "step", "parameters", "n_barcodes", "median_rna_umi", "rna_umi_iqr",
                  "median_genes", "genes_iqr", "median_atac", "atac_iqr", "atac_measure", "note")

qc_summary_file <- function(sample_id, annotation) {
  dir.create("qc_summary", showWarnings = FALSE)
  file.path("qc_summary", paste0(sample_id, ".", annotation, ".qc_steps.tsv"))
}

# one row; umi / genes / atac are per-nucleus vectors for the barcodes left after the step
qc_row <- function(stage, step, parameters, umi, genes, atac = NULL, atac_measure = "", note = "") {
  q <- function(x) {
    x <- as.numeric(x[!is.na(x)])
    if (length(x) == 0) return(c(NA, NA, NA))
    round(as.numeric(stats::quantile(x, c(0.25, 0.5, 0.75))))
  }
  iqr <- function(v) if (is.na(v[2])) NA else paste0(v[1], "-", v[3])
  u <- q(umi); g <- q(genes); a <- if (is.null(atac)) c(NA, NA, NA) else q(atac)
  data.frame(stage = as.character(stage), step = step, parameters = parameters, n_barcodes = length(umi),
             median_rna_umi = u[2], rna_umi_iqr = iqr(u), median_genes = g[2], genes_iqr = iqr(g),
             median_atac = a[2], atac_iqr = iqr(a), atac_measure = if (is.null(atac)) NA else atac_measure,
             note = note, stringsAsFactors = FALSE)
}

# write the rows; with replace_stage, rows of that stage already in the file are replaced and the others kept
qc_write <- function(rows, file, replace_stage = NULL) {
  tab <- do.call(rbind, rows)
  tab[] <- lapply(tab, as.character)
  if (!is.null(replace_stage) && file.exists(file)) {
    old <- read.delim(file, colClasses = "character", check.names = FALSE, na.strings = "")
    old <- old[old$stage != replace_stage, QC_BASE_COLS, drop = FALSE]
    tab <- rbind(old, tab[, QC_BASE_COLS])
  }
  n <- as.numeric(tab$n_barcodes)
  prev <- c(NA, head(n, -1))
  ref <- match("EmptyDrops nuclei calls", tab$step)
  tab$n_removed <- prev - n
  tab$pct_retained_vs_previous <- round(100 * n / prev, 1)
  tab$pct_of_emptydrops_calls <- if (is.na(ref)) NA else ifelse(seq_along(n) >= ref, round(100 * n / n[ref], 1), NA)
  out <- tab[, c("stage", "step", "parameters", "n_barcodes", "n_removed", "pct_retained_vs_previous",
                 "pct_of_emptydrops_calls", setdiff(QC_BASE_COLS, c("stage", "step", "parameters", "n_barcodes")))]
  write.table(out, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  message("QC summary written: ", file)
  invisible(out)
}

# the stage 5 rows, shared by seurat_soupx.R and qc_steps_stage5.R (reconstruction from existing outputs)
#   cells: named list of barcode vectors after each step (barhop, emptydrops, mito, atac_min, atac_high,
#          doublets, umi_cap); umi / genes / atac: named per-barcode vectors; p: parameters list;
#   soup_umi / soup_genes: per-nucleus RNA metrics after SoupX (named by barcode)
qc_stage5_rows <- function(cells, umi, genes, atac, p, soup_umi, soup_genes) {
  am <- "ATAC counts in cellranger-arc peaks"
  r <- function(step, key, params, note = "")
    qc_row(5, step, params, umi[cells[[key]]], genes[cells[[key]]], atac[cells[[key]]], am, note)
  none <- function(x) if (length(x)) paste(x, collapse = ",") else "none"
  soup_cells <- names(soup_umi)
  removed <- umi[soup_cells] - soup_umi
  list(
    r("Barcodes above barhop", "barhop",
      paste0("RNA UMIs >= barhop = ", p$barhop),
      "candidate barcodes, not nuclei; those below 'ambient' form the soup"),
    r("EmptyDrops nuclei calls", "emptydrops",
      paste0("lower = ambient = ", p$ambient, " UMIs; FDR < 0.001 (barcodes < barhop removed first)")),
    r("Mitochondrial filter", "mito",
      paste0("% mito UMIs (all genes) < MT_PCT_MAX = ", p$mt_max, "%"),
      paste0("data-driven cut-off (median + 3 MAD) would be ", round(p$mt_mad, 2), "%; not applied")),
    r("ATAC minimum counts", "atac_min",
      paste0("ATAC counts in peaks > MIN_ATAC_FRAGS = ", p$min_frags)),
    r("ATAC high-count outliers", "atac_high",
      paste0("ATAC counts in peaks <= median + 5 MAD of the nuclei = ", round(p$atac_hi))),
    r("Doublet removal", "doublets",
      paste0("scDblFinder class = singlet; automatic score threshold = ", signif(p$dbl_thr, 3),
             "; expected doublet rate = scDblFinder default"),
      paste0(p$n_dbl, " doublets called (", round(100 * p$n_dbl / length(cells$atac_high), 1), "% of input)",
             if (!is.null(p$dbl_note)) paste0("; ", p$dbl_note) else "")),
    r("High RNA UMI cap", "umi_cap",
      paste0("RNA UMIs < MAX_RNA_UMI = ", p$max_umi)),
    qc_row(5, "SoupX ambient RNA correction",
           paste0("contamination rho = ", round(p$rho, 3), "; soupRange = ", p$barhop, "-", p$ambient,
                  " UMIs; soupQuantile = ", p$soup_quantile, "; germline clusters excluded from the estimate: ",
                  none(p$germline_clusters)),
           soup_umi, soup_genes, atac[soup_cells], am,
           paste0("RNA metrics after correction; median ", round(stats::median(removed)), " UMIs (",
                  round(100 * stats::median(removed / pmax(umi[soup_cells], 1)), 1), "%) removed per nucleus"))
  )
}
