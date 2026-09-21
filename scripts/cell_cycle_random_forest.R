# Early vs late cell-cycle classification with a random forest (adapted from cell_cycle_random_forest.R)
# changes: input/output paths as arguments, output directory created before writing, accuracy computed on the
#          test set actually built, unused 200-cell sub-sampling removed
library(Seurat)
library(Signac)
library(Hmisc)
library(caret)
library(randomForest)

args = commandArgs(trailingOnly=TRUE)
combined_seurat = readRDS(args[1])
out_rds = args[2]
dir.create("cell_cycle_prediction", recursive = T, showWarnings = F)

# training: early and late cell-cycle cells of AB/E/MS/C lineages
train_early = c("earlyABxxx", "earlyEMS", "earlyCx", "earlyMSxx")
train_late = c("ABala", "ABalp", "ABpxa", "EMS", "Ca", "Cp", "MSxa", "MSxp")
test_early = c("earlyABxxxx", "earlyEandMS", "earlyCxx", "earlyMSx")
test_late = c("ABaraa", "ABarpa", "ABpxap", "E", "MS", "Cxa", "Cxp", "MSa", "MSp")

scale_data = combined_seurat@assays$SCT@scale.data
build_set = function(early_types, late_types) {
  early_bc = colnames(combined_seurat)[combined_seurat$cell_type %in% early_types]
  late_bc = colnames(combined_seurat)[combined_seurat$cell_type %in% late_types]
  s = as.data.frame(t(scale_data[, c(early_bc, late_bc), drop = FALSE]))
  s$early_late = c(rep("early", length(early_bc)), rep("late", length(late_bc)))
  s
}
train_set_whole = build_set(train_early, train_late)
test_set_whole = build_set(test_early, test_late)
message("training cells: ", paste(names(table(train_set_whole$early_late)), table(train_set_whole$early_late), collapse = " "))
if (length(unique(train_set_whole$early_late)) < 2) stop("training set needs early and late cells; check cell type annotation")

set.seed(6)
train_rf_all <- train(early_late ~ .,
                      method = "rf",
                      tuneGrid = data.frame(mtry = 50), ntree=100,
                      data = train_set_whole)

if (nrow(test_set_whole) > 0) {
  rf_preds <- predict(train_rf_all, test_set_whole)
  accuracy = mean(rf_preds == test_set_whole$early_late)
  message("test accuracy: ", round(accuracy, 3))
  writeLines(paste("test_accuracy", accuracy), "cell_cycle_prediction/cell_cycle_RF.accuracy.txt")
}

rf_preds_all <- predict(train_rf_all, as.data.frame(t(scale_data)))
combined_seurat$cell_cycle_rf = as.character(rf_preds_all)

write.table(data.frame(barcode = colnames(combined_seurat), cell_cycle = combined_seurat$cell_cycle_rf),
            file = "cell_cycle_prediction/cell_cycle_RF.txt", quote = F, sep = "\t", row.names = F, col.names = F)

pdf("cell_cycle_prediction/cell_cycle_RF.pdf", width = 6, height = 6, useDingbats = F)
print(DimPlot(combined_seurat, reduction = "umapWNN", group.by = "cell_cycle_rf", order = T))
dev.off()

saveRDS(combined_seurat, file = out_rds)
