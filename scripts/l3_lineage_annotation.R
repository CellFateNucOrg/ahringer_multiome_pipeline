# Marker-based tissue/cell-type annotation for L3 C. elegans snMultiome data
# Generic replacement for lineage_specific_annotation.R (no AB/E_MS_C_D embryo splits)
#
# Usage: Rscript scripts/l3_lineage_annotation.R \
#            <seurat_object.rds> <marker_list.txt> <annotation> [taxonomy.tsv]
#
# marker_list.txt  — tab-separated, one row per marker table:
#                    <path_to_marker_tsv>  [<short_label>]
#                    Lines starting with '#' are ignored.
# taxonomy.tsv     — optional; maps raw column names to tissue_type for per-study files
#                    (columns: raw_name, tissue_type, tissue_subtype, canonical_cell_type, canonical_col)
#
# For each marker table the script:
#   1. Scores every WNN cluster by binarised average RNA expression vs. the marker matrix
#   2. Assigns the best-matching cell type (and tissue type) to each cluster
#   3. Saves a per-cell assignment table in cell_type_annotation_all/
#   4. Writes PDFs: binarisation traces, per-table UMAP, combined comparison UMAP
#
# Tissue-type colours are consistent across all UMAPs.

library(ggplot2)
library(dplyr)
library(Seurat)
library(Signac)
library(Hmisc)
library(stringr)
library(multimode)

source("scripts/second_derivative_binarisation.R")
source("scripts/pipeline_utils.R")

plot_dir <- make_plot_dir("l3_lineage_annotation")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: l3_lineage_annotation.R <seurat.rds> <marker_list.txt> <annotation> [taxonomy.tsv]")
seurat_object_file <- args[1]
marker_list_file   <- args[2]
annotation         <- args[3]
taxonomy_file      <- if (length(args) >= 4) args[4] else ""

# ── Consistent tissue-type colour palette ─────────────────────────────────────
TISSUE_COLOURS <- c(
  "Neurons"      = "#5B5EA6",  # indigo-blue
  "Muscle"       = "#E63946",  # red
  "Intestine"    = "#57CC99",  # green
  "Hypodermis"   = "#F4D35E",  # yellow
  "Pharynx"      = "#F77F00",  # orange
  "Germline"     = "#C77DFF",  # violet
  "Glia"         = "#48CAE4",  # sky blue
  "Coelomocytes" = "#1B998B",  # teal
  "Excretory"    = "#E9C46A",  # sand
  "Gonad"        = "#FFA8A8",  # light pink
  "unassigned"   = "#AAAAAA",  # grey
  "Other"        = "#CCCCCC"   # light grey
)

# ── Scoring function (from lineage_specific_annotation.R) ─────────────────────
scoring_markers <- function(marker_col, cluster_col) {
  # marker_col, cluster_col : numeric vectors of same length, values in {0,1,2}
  score  <- length(marker_col)
  malus  <- score * 0.1
  kill   <- malus  * 5
  for (i in seq_along(marker_col)) {
    d <- abs(marker_col[i] - cluster_col[i])
    if      (d == 2) score <- score - kill
    else if (d == 1) score <- score - malus
  }
  score
}

# ── Binarisation with fallback ─────────────────────────────────────────────────
binarise_gene <- function(values, gene_name) {
  tryCatch({
    loc <- locmodes(values, mod0 = 2, display = FALSE)
    second_derivative_binarization(values, loc$cbw$bw, loc$locations, gene_name)
  }, error = function(e) {
    warning("binarisation failed for ", gene_name, " (", conditionMessage(e),
            "); using 25%/60% of max as fallback thresholds")
    m <- max(values)
    if (m <= 0) c(Inf, Inf) else c(0.25 * m, 0.6 * m)
  })
}

# ── Tissue-type extraction ─────────────────────────────────────────────────────
# Canonical columns:  "Neurons.Dopaminergic.CEP" -> "Neurons"
# Raw columns:        look up in taxonomy table
extract_tissue_type <- function(col_name, taxonomy = NULL) {
  if (grepl("\\.", col_name)) {
    return(strsplit(col_name, "\\.")[[1]][1])
  }
  if (!is.null(taxonomy) && col_name %in% taxonomy$raw_name) {
    return(taxonomy$tissue_type[match(col_name, taxonomy$raw_name)])
  }
  "Other"
}

# ── Load taxonomy (optional) ───────────────────────────────────────────────────
taxonomy <- NULL
if (taxonomy_file != "" && file.exists(taxonomy_file)) {
  taxonomy <- read.table(taxonomy_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  message("Taxonomy loaded: ", nrow(taxonomy), " entries")
}

# ── Read marker list ───────────────────────────────────────────────────────────
raw_lines <- readLines(marker_list_file)
raw_lines <- raw_lines[!grepl("^\\s*#", raw_lines) & nzchar(trimws(raw_lines))]
marker_entries <- do.call(rbind, lapply(raw_lines, function(ln) {
  parts <- strsplit(ln, "\t")[[1]]
  path  <- trimws(parts[1])
  label <- if (length(parts) >= 2) trimws(parts[2]) else ""
  if (label == "") label <- tools::file_path_sans_ext(basename(path))
  data.frame(path = path, label = label, stringsAsFactors = FALSE)
}))
message("Marker list: ", nrow(marker_entries), " tables")

# ── Load Seurat object ─────────────────────────────────────────────────────────
message("Loading Seurat object: ", seurat_object_file)
combined_seurat <- readRDS(seurat_object_file)
n_clusters <- length(levels(combined_seurat$seurat_clusters))
message("Cells: ", ncol(combined_seurat), "  WNN clusters: ", n_clusters)

# ── Process each marker table ──────────────────────────────────────────────────
all_tissue_meta <- list()   # for summary plot
out_dir <- "cell_type_annotation_all"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (row_i in seq_len(nrow(marker_entries))) {
  mpath  <- marker_entries$path[row_i]
  mlabel <- marker_entries$label[row_i]
  safe_label <- gsub("[^A-Za-z0-9._-]", "_", mlabel)
  message("\n=== ", mlabel, " ===")

  if (!file.exists(mpath)) {
    warning("Marker file not found, skipping: ", mpath)
    next
  }

  markers_table <- read.table(mpath, header = TRUE, sep = "\t",
                               row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  message("  Table: ", nrow(markers_table), " genes x ", ncol(markers_table), " cell types")

  rna_genes     <- rownames(combined_seurat[["RNA"]])
  markers_table <- markers_table[rownames(markers_table) %in% rna_genes, , drop = FALSE]
  if (nrow(markers_table) == 0) {
    warning("No marker genes found in RNA assay for '", mlabel, "', skipping")
    next
  }
  message("  Genes matched to RNA assay: ", nrow(markers_table))

  # Average expression per WNN cluster
  Idents(combined_seurat) <- "seurat_clusters"
  avg_exp     <- as.data.frame(AverageExpression(combined_seurat, slot = "data",
                                                  assays = "RNA", features = rownames(markers_table)))
  avg_exp     <- avg_exp[rownames(markers_table), , drop = FALSE]
  avg_exp_log <- log10(avg_exp + 1)

  # Column names: "RNA.0", "RNA.1", ... -> strip prefix to get cluster IDs
  cluster_ids <- gsub("^RNA\\.", "", colnames(avg_exp_log))

  # Binarise expression per gene
  pdf(file = file.path(plot_dir, paste0(safe_label, ".binarisation.", annotation, ".pdf")),
      width = 8, height = 8)
  avg_exp_quant <- avg_exp_log
  for (gene_i in seq_len(nrow(avg_exp_log))) {
    values     <- as.numeric(avg_exp_log[gene_i, ])
    thresholds <- binarise_gene(values, rownames(avg_exp_log)[gene_i])
    avg_exp_quant[gene_i, values <  thresholds[1]]                           <- 0
    avg_exp_quant[gene_i, values >= thresholds[1] & values < thresholds[2]] <- 1
    avg_exp_quant[gene_i, values >= thresholds[2]]                           <- 2
  }
  dev.off()

  # Score clusters vs. cell types
  n_ct      <- ncol(markers_table)
  n_cl      <- length(cluster_ids)
  n_markers <- nrow(markers_table)
  scoring   <- matrix(NA_real_, nrow = n_ct, ncol = n_cl,
                      dimnames = list(colnames(markers_table), cluster_ids))
  for (cl_i in seq_len(n_cl)) {
    cl_quant <- as.numeric(avg_exp_quant[, cl_i])
    for (ct_i in seq_len(n_ct)) {
      scoring[ct_i, cl_i] <- scoring_markers(markers_table[, ct_i], cl_quant)
    }
  }

  write.table(scoring,
              file = file.path(out_dir, paste0(safe_label, ".cluster_scoring.", annotation, ".txt")),
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

  # Assign clusters
  cluster_assignment  <- setNames(rep(NA_character_, n_cl), cluster_ids)
  cluster_tissue_type <- setNames(rep("unassigned",  n_cl), cluster_ids)

  for (cl_i in seq_len(n_cl)) {
    scores     <- scoring[, cl_i]
    best_score <- max(scores)
    best_idx   <- which(scores == best_score)
    scores_ord <- sort(scores, decreasing = TRUE)
    second_best <- if (length(scores_ord) >= 2) scores_ord[2] else -Inf

    assigned <- NA_character_
    if (best_score >= n_markers * 0.7 && length(best_idx) == 1) {
      assigned <- names(best_idx)
    } else if (best_score > 0 && (best_score - second_best) > n_markers * 0.5) {
      assigned <- names(best_idx[1])
    }
    cluster_assignment [cluster_ids[cl_i]] <- assigned
    if (!is.na(assigned))
      cluster_tissue_type[cluster_ids[cl_i]] <- extract_tissue_type(assigned, taxonomy)
  }

  assigned_n <- sum(!is.na(cluster_assignment))
  message("  Clusters assigned: ", assigned_n, " / ", n_cl)
  if (assigned_n < n_cl / 2)
    warning("  Only ", assigned_n, " / ", n_cl, " clusters assigned — consider using a larger/more specific marker table")

  # Map to cells
  cell_cl_id       <- as.character(combined_seurat$seurat_clusters)
  cell_assignment  <- cluster_assignment [cell_cl_id]
  cell_tissue_type <- cluster_tissue_type[cell_cl_id]
  cell_assignment [is.na(cell_assignment)]  <- "unassigned"
  cell_tissue_type[is.na(cell_tissue_type)] <- "unassigned"

  meta_ct <- paste0(safe_label, "_cell_type")
  meta_tt <- paste0(safe_label, "_tissue_type")
  combined_seurat@meta.data[[meta_ct]] <- cell_assignment
  combined_seurat@meta.data[[meta_tt]] <- cell_tissue_type
  all_tissue_meta[[mlabel]] <- list(meta_col = meta_tt, tissue_vec = cell_tissue_type)

  # Save per-cell assignment table
  out_df <- data.frame(cell_type   = cell_assignment,
                       tissue_type = cell_tissue_type,
                       row.names   = colnames(combined_seurat))
  write.table(out_df,
              file = file.path(out_dir, paste0(safe_label, ".cell_assignment.", annotation, ".txt")),
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

  # ── UMAPs ──────────────────────────────────────────────────────────────────
  tt_present <- sort(unique(cell_tissue_type))
  col_map    <- TISSUE_COLOURS[tt_present]
  missing    <- is.na(col_map)
  if (any(missing)) {
    col_map[missing] <- grDevices::colorRampPalette(c("#E8E8E8", "#B0B0B0"))(sum(missing))
    names(col_map)[missing] <- tt_present[missing]
  }

  pdf(file = file.path(plot_dir, paste0(safe_label, ".annotation.UMAP.", annotation, ".pdf")),
      width = 16, height = 7, useDingbats = FALSE)

  p_tissue <- DimPlot(combined_seurat, reduction = "umapWNN", group.by = meta_tt,
                      label = TRUE, label.size = 3, order = TRUE, repel = TRUE) +
    scale_color_manual(values = col_map) +
    ggtitle(paste0(mlabel, " — tissue type")) +
    theme(legend.text  = element_text(size = 8),
          plot.title   = element_text(size = 12))

  p_cell <- DimPlot(combined_seurat, reduction = "umapWNN", group.by = meta_ct,
                    label = TRUE, label.size = 2.5, order = TRUE, repel = TRUE) +
    ggtitle(paste0(mlabel, " — cell type")) +
    NoLegend()

  print(p_tissue | p_cell)
  dev.off()
  message("  Saved UMAP: ", safe_label, ".annotation.UMAP.", annotation, ".pdf")
}

# ── Summary comparison UMAP ────────────────────────────────────────────────────
if (length(all_tissue_meta) >= 2) {
  message("\nGenerating summary comparison UMAP (", length(all_tissue_meta), " tables)...")
  n_tables <- length(all_tissue_meta)
  n_cols   <- ceiling(sqrt(n_tables))
  n_rows   <- ceiling(n_tables / n_cols)

  pdf(file = file.path(plot_dir, paste0("all_tables.tissue_type.UMAP.", annotation, ".pdf")),
      width = 8 * n_cols, height = 7 * n_rows, useDingbats = FALSE)

  plot_list <- lapply(names(all_tissue_meta), function(lbl) {
    info       <- all_tissue_meta[[lbl]]
    tt_present <- sort(unique(info$tissue_vec))
    col_map    <- TISSUE_COLOURS[tt_present]
    missing    <- is.na(col_map)
    if (any(missing)) {
      col_map[missing] <- grDevices::colorRampPalette(c("#E8E8E8","#B0B0B0"))(sum(missing))
      names(col_map)[missing] <- tt_present[missing]
    }
    DimPlot(combined_seurat, reduction = "umapWNN", group.by = info$meta_col,
            label = TRUE, label.size = 3, order = TRUE, repel = TRUE) +
      scale_color_manual(values = col_map) +
      ggtitle(lbl) +
      theme(plot.title    = element_text(size = 10),
            legend.text   = element_text(size = 7),
            legend.key.size = unit(0.4, "cm"))
  })

  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    print(wrap_plots(plot_list, ncol = n_cols))
  } else {
    for (p in plot_list) print(p)
  }
  dev.off()
  message("Saved: all_tables.tissue_type.UMAP.", annotation, ".pdf")
}

message("\nAll done.")
message("  Plots      -> ", plot_dir, "/")
message("  Assignments-> ", out_dir, "/")
