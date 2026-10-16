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

# ── Scoring function ───────────────────────────────────────────────────────────
# Returns a NORMALISED score in (-inf, 1].
# Each gene contributes: +1 (perfect match), -0.5 (off by 1), -2.5 (off by 2).
# Dividing by n keeps the score on the same scale regardless of marker-table size,
# which is critical for large tables (100s–1000s of genes) where the original
# absolute penalty would make every score deeply negative.
scoring_markers <- function(marker_col, cluster_col) {
  n <- length(marker_col)
  if (n == 0) return(0)
  score <- 0
  for (i in seq_along(marker_col)) {
    d <- abs(marker_col[i] - cluster_col[i])
    if      (d == 0) score <- score + 1.0   # perfect
    else if (d == 1) score <- score - 0.5   # one step off
    else if (d == 2) score <- score - 2.5   # two steps off (same 5:1 kill:malus ratio)
  }
  score / n   # range: -2.5 (all kills) to 1.0 (all perfect)
}

# ── Binarisation with fallback ─────────────────────────────────────────────────
# On failure, return c(Inf, Inf) so all values → 0 (conservative: no false-high).
# The previous fallback c(0.25*m, 0.6*m) with tiny m caused spurious value=2
# assignments that killed scores via the penalty function.
binarise_gene <- function(values, gene_name) {
  tryCatch({
    loc <- locmodes(values, mod0 = 2, display = FALSE)
    second_derivative_binarization(values, loc$cbw$bw, loc$locations, gene_name)
  }, error = function(e) {
    warning("binarisation failed for ", gene_name, " (", conditionMessage(e),
            "); treating as unexpressed (Inf thresholds)")
    c(Inf, Inf)   # all values → 0; avoids spurious kills from relative thresholds
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

# ── Marker gene names -> gene names of this annotation (via WBGeneID) ──────────
# The marker tables use current WormBase gene labels, but the count matrix uses the names of the
# WormBase release of this run (config WS, e.g. WS285). Genes renamed in between (e.g. spig-9 is
# F53F4.13 in WS285) would silently drop out, so labels that are not names of this annotation are
# translated label -> WBGeneID (tissueMarkers evidence / consensus tables) -> annotation name.
id_tables <- strsplit(env_chr("MARKER_GENE_ID_TABLES",
  "tissueMarkers/evidence/marker_evidence_long_tissue.csv tissueMarkers/consensus/markers_consensus_min2sources.csv"),
  "[ ,]+")[[1]]
label_to_id <- character(0)
for (f in id_tables[file.exists(id_tables)]) {
  id_tab <- read.csv(f, stringsAsFactors = FALSE)
  if (all(c("gene_label", "WBGeneID") %in% colnames(id_tab)))
    label_to_id <- c(label_to_id, setNames(id_tab$WBGeneID, id_tab$gene_label))
}
label_to_id <- label_to_id[!duplicated(names(label_to_id))]
annot_gene_names <- read.table(paste0(canon_prefix(), ".gene_name.txt"), col.names = c("wbgene", "gene_name"),
                               stringsAsFactors = FALSE)
id_to_name <- setNames(annot_gene_names$gene_name, annot_gene_names$wbgene)
message("Gene-name translation: ", length(label_to_id), " marker labels with WBGeneIDs, ",
        env_chr("WS", "WS285"), " annotation with ", nrow(annot_gene_names), " genes")
to_annotation_names <- function(g) {
  need <- which(!(g %in% annot_gene_names$gene_name) & g %in% names(label_to_id))
  new  <- unname(id_to_name[label_to_id[g[need]]])
  ok   <- !is.na(new)
  g[need[ok]] <- new[ok]
  g
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

  # translate marker labels to the gene names of this annotation (see to_annotation_names above)
  orig_names <- rownames(markers_table)
  new_names  <- to_annotation_names(orig_names)
  renamed    <- which(new_names != orig_names)
  if (length(renamed) > 0) {
    message("  Renamed to ", env_chr("WS", "WS285"), " names via WBGeneID: ", length(renamed), " (",
            paste(head(paste0(orig_names[renamed], "->", new_names[renamed]), 8), collapse = ", "),
            if (length(renamed) > 8) ", ..." else "", ")")
  }
  keep <- !duplicated(new_names)
  if (any(!keep)) message("  Dropped ", sum(!keep), " duplicate gene(s) after renaming")
  markers_table <- markers_table[keep, , drop = FALSE]
  rownames(markers_table) <- new_names[keep]
  # names shown in plots: the common name wherever one exists (marker label if it is a common name,
  # else the annotation name), e.g. Y75B7AL.2 is shown as amn-1 and F57F5.1 as cpr-9
  display_name <- setNames(ifelse(is_common_name(orig_names[keep]) | !is_common_name(new_names[keep]),
                                  orig_names[keep], new_names[keep]),
                           new_names[keep])

  rna_genes     <- rownames(combined_seurat[["RNA"]])
  markers_table <- markers_table[rownames(markers_table) %in% rna_genes, , drop = FALSE]
  if (nrow(markers_table) == 0) {
    warning("No marker genes found in RNA assay for '", mlabel, "', skipping")
    next
  }
  message("  Genes matched to RNA assay: ", nrow(markers_table))

  # Average expression per WNN cluster
  Idents(combined_seurat) <- "seurat_clusters"
  avg_exp     <- as.data.frame(AverageExpression(combined_seurat, layer = "data",
                                                  assays = "RNA", features = rownames(markers_table)))
  avg_exp     <- avg_exp[rownames(markers_table), , drop = FALSE]
  avg_exp_log <- log10(avg_exp + 1)

  # Column names: "RNA.0", "RNA.1", ... -> strip prefix to get cluster IDs
  cluster_ids <- gsub("^RNA\\.g?", "", colnames(avg_exp_log))

  # Binarise expression per gene
  pdf(file = file.path(plot_dir, paste0(safe_label, ".binarisation.", annotation, ".pdf")),
      width = 8, height = 8)
  avg_exp_quant <- avg_exp_log
  for (gene_i in seq_len(nrow(avg_exp_log))) {
    values     <- as.numeric(avg_exp_log[gene_i, ])
    thresholds <- binarise_gene(values, display_name[[rownames(avg_exp_log)[gene_i]]])
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
    # Thresholds on the NORMALISED scale (scoring_markers now returns score/n_genes).
    # Condition 1: best score ≥ 0.5 AND unique winner  → high-confidence assignment
    # Condition 2: best score > 0   AND gap > 0.10     → weaker but distinct winner
    # These replace the old absolute thresholds (n_markers*0.7, n_markers*0.5) that
    # were only valid for small (10–20 gene) marker tables.
    if (best_score >= 0.5 && length(best_idx) == 1) {
      assigned <- names(best_idx)
    } else if (best_score > 0.0 && (best_score - second_best) > 0.10) {
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
