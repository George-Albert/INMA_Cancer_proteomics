suppressPackageStartupMessages({
  library(limma)
  library(preprocessCore)
})

source("src/utils.R")
source("src/statistics.R")

# Step 4: run limma differential expression across mean-filter thresholds.
input_dir <- "Inputs"
output_dir <- "Outputs"

metadata <- read.table(file.path(input_dir, "txt", "metadata_filtered.txt"))
reads_spec <- read.table(file.path(input_dir, "txt", "reads_spec_filtered.txt"))

th_mean_vec <- c(0.2, 0.5, 1)
p_value_threshold <- 0.05
rm_ext <- "both"

design <- build_design(metadata)
contrast_matrix <- build_default_contrasts(design)

create_dir(file.path(output_dir, "plot_SA"))
create_dir(file.path(input_dir, "003_DEG"))
create_dir(file.path(input_dir, "004_Significative_DEG"))

for (th_mean in th_mean_vec) {
  exp_norm <- normalize_expression(reads_spec, rm_ext = rm_ext)
  means <- apply(exp_norm, 1, mean)
  reads_filtered <- reads_spec[which(means > th_mean), , drop = FALSE]

  exp_norm_threshold <- normalize_expression(reads_filtered, rm_ext = rm_ext)
  fit <- limma::lmFit(exp_norm_threshold, design)
  fit_bayes <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

  tag <- build_threshold_tag(th_mean)
  sa_plot_name <- paste0("Plot_SA_", tag, ".pdf")
  grDevices::pdf(file.path(output_dir, "plot_SA", sa_plot_name), width = 8, height = 6)
  limma::plotSA(fit_bayes, xlab = "Average log-expression", ylab = "sqrt(sigma)", pch = 19, cex = 0.6)
  grDevices::dev.off()

  total_fit <- limma::contrasts.fit(fit, contrast_matrix)
  total_fit <- limma::eBayes(total_fit)
  n_genes <- nrow(reads_filtered)

  for (contrast_name in colnames(contrast_matrix)) {
    all_deg <- limma::topTable(total_fit, coef = contrast_name, number = n_genes, adjust.method = "BH")
    sig_deg <- limma::topTable(total_fit, coef = contrast_name, number = n_genes, adjust.method = "BH", p.value = p_value_threshold)

    write.table(all_deg, file = file.path(input_dir, "003_DEG", paste0("deg_", contrast_name, tag, ".txt")))
    write.table(sig_deg, file = file.path(input_dir, "004_Significative_DEG", paste0("deg_", contrast_name, tag, ".txt")))
  }
}
