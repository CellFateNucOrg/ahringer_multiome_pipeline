# Final cell type assignment for L3 C. elegans snMultiome data
# L3-adapted replacement for the embryo final_cell_type_assignment.R
#
# Changes vs. embryo version:
#   * No AB / E_MS_C_D lineage tables — reads l3_lineage_annotation output files instead
#   * No cpg-3 two-cell, prg-2 Z2/Z3, P-granule P-cell, or m1m2 logic
#   * No embryonic AB/MS/C/D lineage grouping
#   * "cell_lineage" is set to tissue_type (broadest grouping available at L3)
#   * Primary annotation chosen via PRIMARY_ANNOTATION_LABEL env var or command-line arg
#
# Usage: Rscript scripts/final_cell_type_assignment.R \
#            <seurat_object.rds> <annotation> [primary_marker_label]
#
# primary_marker_label: short label used by l3_lineage_annotation (e.g. "shared_min2studies")
#                       Matches the second column of cellMarkers/marker_list.txt.
#                       Defaults to PRIMARY_ANNOTATION_LABEL env var, then to the first
#                       available label found in cell_type_annotation_all/.

library(ggplot2)
library(dplyr)
library(Seurat)
library(Signac)

source("scripts/pipeline_utils.R")

plot_dir <- make_plot_dir("final_cell_types")

args       <- commandArgs(trailingOnly = TRUE)
seurat_rds <- args[1]
annotation <- args[2]
primary_label_arg <- if (length(args) >= 3) args[3] else ""

# ── Tissue-type colour palette (matches l3_lineage_annotation.R) ───────────────
TISSUE_COLOURS <- c(
  "Neurons"      = "#5B5EA6",
  "Muscle"       = "#E63946",
  "Intestine"    = "#57CC99",
  "Hypodermis"   = "#F4D35E",
  "Pharynx"      = "#F77F00",
  "Germline"     = "#C77DFF",
  "Glia"         = "#48CAE4",
  "Coelomocytes" = "#1B998B",
  "Excretory"    = "#E9C46A",
  "Gonad"        = "#FFA8A8",
  "unassigned"   = "#AAAAAA",
  "Other"        = "#CCCCCC"
)

# ── Helpers ────────────────────────────────────────────────────────────────────
pick_primary_label <- function(annotation, arg_label) {
  candidates <- c(arg_label, Sys.getenv("PRIMARY_ANNOTATION_LABEL", ""))
  for (lbl in candidates) {
    if (nzchar(lbl)) {
      f <- file.path("cell_type_annotation_all",
                     paste0(lbl, ".cell_assignment.", annotation, ".txt"))
      if (file.exists(f)) {
        message("Primary annotation label: '", lbl, "' (explicit)")
        return(lbl)
      }
      warning("Specified primary label '", lbl, "' has no matching file, trying auto-detect")
    }
  }
  pattern <- paste0("\\.cell_assignment\\.", annotation, "\\.txt$")
  files   <- list.files("cell_type_annotation_all", pattern = pattern, full.names = FALSE)
  if (length(files) == 0) stop("No cell_assignment files found in cell_type_annotation_all/")
  lbl <- sub(paste0("\\.cell_assignment\\.", annotation, "\\.txt$"), "", files[1])
  message("Primary annotation label: '", lbl, "' (auto-detected)")
  lbl
}

colour_umap <- function(obj, group_by, colours, title, label_size = 3) {
  vals    <- sort(unique(as.character(obj@meta.data[[group_by]])))
  col_map <- colours[vals]
  missing <- is.na(col_map)
  if (any(missing)) {
    col_map[missing] <- grDevices::colorRampPalette(c("#E8E8E8","#B0B0B0"))(sum(missing))
    names(col_map)[missing] <- vals[missing]
  }
  DimPlot(obj, reduction = "umapWNN", group.by = group_by,
          label = TRUE, label.size = label_size, order = TRUE, repel = TRUE) +
    scale_color_manual(values = col_map) +
    ggtitle(title) +
    theme(legend.text  = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"),
          plot.title   = element_text(size = 11))
}

# ── Load Seurat object ─────────────────────────────────────────────────────────
message("Loading Seurat object: ", seurat_rds)
combined_seurat <- readRDS(seurat_rds)
message("Cells: ", ncol(combined_seurat))

out_prefix <- paste0("seurat_objects/wt_all.all_samples.all_cells.", annotation)

# ── Choose primary annotation label ────────────────────────────────────────────
primary_label <- pick_primary_label(annotation, primary_label_arg)
safe_label    <- gsub("[^A-Za-z0-9._-]", "_", primary_label)

ct_file <- file.path("cell_type_annotation_all",
                     paste0(safe_label, ".cell_assignment.", annotation, ".txt"))
message("Reading primary cell assignment: ", ct_file)
primary_assign <- read.table(ct_file, header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE, row.names = 1)

# ── Assign cell_type and cell_lineage ──────────────────────────────────────────
# cell_type    = fine-grained assignment (e.g. Neurons.Interneuron.AIY or "unassigned")
# cell_lineage = tissue-level grouping (e.g. Neurons, Muscle) -- replaces embryo lineages at L3
matched_idx  <- match(colnames(combined_seurat), rownames(primary_assign))
cell_type    <- ifelse(is.na(matched_idx), "unassigned", primary_assign$cell_type[matched_idx])
cell_lineage <- ifelse(is.na(matched_idx), "unassigned", primary_assign$tissue_type[matched_idx])
cell_type   [is.na(cell_type)]    <- "unassigned"
cell_lineage[is.na(cell_lineage)] <- "unassigned"

combined_seurat$cell_type    <- cell_type
combined_seurat$cell_lineage <- cell_lineage

n_assigned <- sum(cell_type != "unassigned")
message("Assigned: ", n_assigned, " / ", ncol(combined_seurat), " cells")
message("Tissue types: ", paste(sort(unique(cell_lineage)), collapse = ", "))
print(table(cell_lineage))

# ── Apply UNASSIGN_CELL_TYPES (set in config.sh or env) ───────────────────────
unassign_types <- strsplit(env_chr("UNASSIGN_CELL_TYPES", ""), "[ ,]+")[[1]]
unassign_types <- unassign_types[nzchar(unassign_types)]
if (length(unassign_types) > 0) {
  message("Unassigning: ", paste(unassign_types, collapse = ", "))
  combined_seurat$cell_type   [combined_seurat$cell_type    %in% unassign_types] <- "unassigned"
  combined_seurat$cell_lineage[combined_seurat$cell_lineage %in% unassign_types] <- "unassigned"
}

# ── Merge all other annotation labels into metadata (for comparison) ───────────
all_files <- list.files("cell_type_annotation_all",
                        pattern = paste0("\\.cell_assignment\\.", annotation, "\\.txt$"),
                        full.names = FALSE)
for (fn in all_files) {
  lbl  <- sub(paste0("\\.cell_assignment\\.", annotation, "\\.txt$"), "", fn)
  slbl <- gsub("[^A-Za-z0-9._-]", "_", lbl)
  if (slbl == safe_label) next
  tbl  <- tryCatch(
    read.table(file.path("cell_type_annotation_all", fn), header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, row.names = 1),
    error = function(e) { warning("Could not read ", fn, ": ", conditionMessage(e)); NULL }
  )
  if (is.null(tbl)) next
  idx <- match(colnames(combined_seurat), rownames(tbl))
  combined_seurat@meta.data[[paste0(slbl, "_cell_type")]]   <- ifelse(is.na(idx), "unassigned", tbl$cell_type[idx])
  combined_seurat@meta.data[[paste0(slbl, "_tissue_type")]] <- ifelse(is.na(idx), "unassigned", tbl$tissue_type[idx])
}

# ── UMAPs ──────────────────────────────────────────────────────────────────────
message("Writing UMAPs...")
pdf(file = file.path(plot_dir, paste0("wt.all_samples.all_cells.", annotation, ".UMAP.pdf")),
    width = 16, height = 8, useDingbats = FALSE)
p_lin <- colour_umap(combined_seurat, "cell_lineage", TISSUE_COLOURS,
                     paste0(primary_label, " — tissue type"), label_size = 3.5)
p_cl  <- colour_umap(combined_seurat, "cell_type",    TISSUE_COLOURS,
                     paste0(primary_label, " — cell type"),   label_size = 2.5)
print(p_lin | p_cl)
dev.off()

# ── Fragment BED files per group (needed by stage 9 peaks_idr) ─────────────────
message("Splitting fragments by cell_lineage and cell_type...")
DefaultAssay(combined_seurat) <- "WNN"
for (cell_annotation in c("cell_lineage", "cell_type")) {
  d <- file.path("cluster_specific_peaks", "round_2_all_samples", annotation, cell_annotation, "bed")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(file.path(d, "*.bed"))
  SplitFragments(combined_seurat, group.by = cell_annotation, outdir = d, append = TRUE)
  write_cells_per_group(combined_seurat, cell_annotation,
                        paste0(out_prefix, ".cells_per_", cell_annotation, ".txt"))
  message("  ", cell_annotation, " -> ", d)
}

# ── Save RDS and table ─────────────────────────────────────────────────────────
message("Saving: ", out_prefix, ".rds")
saveRDS(combined_seurat, file = paste0(out_prefix, ".rds"))
write.table(combined_seurat[[]], file = paste0(out_prefix, ".table.txt"),
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

message("\nAll done.")
message("  RDS   -> ", out_prefix, ".rds")
message("  Table -> ", out_prefix, ".table.txt")
message("  Plots -> ", plot_dir, "/")
print(table(combined_seurat$cell_type))
