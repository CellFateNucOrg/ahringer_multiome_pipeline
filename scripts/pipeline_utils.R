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

key_lineage_markers <- c("med-1", "tbx-40", "tbx-35", "tbx-33", "end-3", "hnd-1", "pes-10", "vet-2", "ceh-51",
                         "oma-1", "oma-2", "lsl-1", "skr-7", "ccch-2", "ref-1", "cpg-3")
feature_colours <- c("grey90", "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred")

plot_key_markers <- function(obj, reduction, file) {
  pdf(file = file, width = 8, height = 8)
  for (g in intersect(key_lineage_markers, rownames(obj))) {
    print(FeaturePlot(obj, order = TRUE, features = g, reduction = reduction) & scale_colour_gradientn(colours = feature_colours))
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
