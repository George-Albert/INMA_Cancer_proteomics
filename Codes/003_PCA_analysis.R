suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(preprocessCore)
})

source("src/utils.R")
source("src/statistics.R")
source("src/visualization.R")

# Step 3: run PCA before and after removing known outlier samples.
input_dir <- "Inputs"
output_dir <- "Outputs"
pca_dir <- "PCA"

feature_data <- read.table(file.path(input_dir, "txt", "feature_data.txt"))
metadata <- read.table(file.path(input_dir, "txt", "metadata.txt"))
reads_spec <- read.table(file.path(input_dir, "txt", "reads_spec.txt"))

if ("PEG_4h_D238_8" %in% colnames(reads_spec)) {
  reads_spec[, "PEG_4h_D238_8"] <- NULL
  metadata <- metadata[rownames(metadata) != "PEG_4h_D238_8", ]
}

reads_spec[is.na(reads_spec)] <- 0
exp_norm <- normalize_expression(reads_spec, rm_ext = "both")
means <- apply(exp_norm, 1, mean)
reads_filtered <- reads_spec[which(means > 1), , drop = FALSE]

run_and_export_pca <- function(reads, metadata, table_name, plot_name) {
  exp_norm_local <- normalize_expression(reads, rm_ext = "both")
  pca <- prcomp(t(exp_norm_local))
  pca_summary <- summary(pca)$importance

  pca_df <- as.data.frame(pca$x)
  combined_data <- dplyr::bind_cols(metadata, pca_df[, 1:6])

  create_dir(file.path(output_dir, pca_dir, "PCA_tables"))
  write.table(combined_data, file = file.path(output_dir, pca_dir, "PCA_tables", paste0(table_name, ".txt")))

  p <- build_pca_plot(combined_data, pca_summary[2, 1], pca_summary[2, 2])

  create_dir(file.path(output_dir, pca_dir))
  grDevices::pdf(file = file.path(output_dir, pca_dir, paste0(plot_name, ".pdf")), width = 6, height = 5)
  print(p)
  grDevices::dev.off()
}

run_and_export_pca(reads_filtered, metadata, "PCA_table", "PCA_all_particles")

outlier_ids <- c(14, 9, 4)
metadata_wo_outlier <- metadata[!(metadata$Sample.Order %in% outlier_ids), , drop = FALSE]
reads_wo_outlier <- reads_filtered[, colnames(reads_filtered) %in% rownames(metadata_wo_outlier), drop = FALSE]

run_and_export_pca(reads_wo_outlier, metadata_wo_outlier, "PCA_table_wo_outlier", "PCA_all_particles_wo_outlier")

write.table(metadata_wo_outlier, file.path(input_dir, "txt", "metadata_filtered.txt"))
write.table(reads_wo_outlier, file.path(input_dir, "txt", "reads_spec_filtered.txt"))
