# Helper (new): barcode -> cell type tables per sample (for splitting STARsolo BAMs, replaces the
# 'cut -f 1,14 ... | grep sample' commands of the paper) and stage -> cell type table for stage ATAC tracks
suppressPackageStartupMessages(library(Seurat))
source("scripts/pipeline_utils.R")

args = commandArgs(trailingOnly = TRUE)
obj = readRDS(args[1])
out_dir = args[2]
dir.create(file.path(out_dir, "rna_barcodes"), recursive = TRUE, showWarnings = FALSE)

md = obj@meta.data
for (id in unique(md$exp)) {
  sel = md$exp == id
  bc = sub(paste0("^", id, "_"), "", rownames(md)[sel])
  write.table(data.frame(bc, md$cell_type[sel]), file = file.path(out_dir, "rna_barcodes", paste0(id, ".barcodes.txt")),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}
stage_tab = unique(data.frame(stage = md$stage, cell_type = md$cell_type))
write.table(stage_tab[order(stage_tab$stage, stage_tab$cell_type), ], file = file.path(out_dir, "stage_cell_types.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
