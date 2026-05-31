suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
})

source("src/utils.R")
source("src/io.R")
source("src/preprocessing.R")

# Step 2: consolidate quality-annotated samples into reads/metadata/feature tables.
paths <- list(
  input = "Inputs",
  quality = "Inputs/002_Quality_annotated",
  txt = "Inputs/txt",
  xlsx = "Inputs/xlsx"
)

create_dir(paths$txt)
create_dir(paths$xlsx)

sample_files <- list_excel_files(paths$quality)
sample_data <- lapply(sample_files, function(x) openxlsx::read.xlsx(file.path(paths$quality, x)))

sample_names <- vapply(sample_files, build_sample_variable_name, character(1))
names(sample_data) <- sample_names

columns_to_keep <- c("Gene.Name", "accession", "description", "spec.count", "EMPAI")
selected <- lapply(seq_along(sample_data), function(i) {
  df <- sample_data[[i]][, columns_to_keep]
  colnames(df)[colnames(df) != "accession"] <- paste0(colnames(df)[colnames(df) != "accession"], "_", names(sample_data)[i])
  df
})

full_data_df <- Reduce(function(x, y) merge(x, y, by = "accession", all = TRUE, sort = FALSE), selected)

gene_name_columns <- grep("^Gene.Name_", colnames(full_data_df), value = TRUE)
description_columns <- grep("^description_", colnames(full_data_df), value = TRUE)
spec_columns <- grep("^spec.count_", colnames(full_data_df), value = TRUE)
empai_columns <- grep("^EMPAI_", colnames(full_data_df), value = TRUE)

full_data_df$Gene.Name <- apply(full_data_df[, gene_name_columns], 1, first_non_na)
full_data_df$Description <- apply(full_data_df[, description_columns], 1, first_non_na)

reads_spec <- full_data_df[, c("Gene.Name", "accession", spec_columns)]
reads_empai <- full_data_df[, c("Gene.Name", "accession", empai_columns)]

colnames(reads_spec) <- sub("^spec.count_", "", colnames(reads_spec))
colnames(reads_empai) <- sub("^EMPAI_", "", colnames(reads_empai))

duplicate_indices <- which(duplicated(reads_spec$Gene.Name))
if (length(duplicate_indices) > 0) {
  reads_spec$Gene.Name[duplicate_indices] <- reads_spec$accession[duplicate_indices]
}

rownames(reads_spec) <- reads_spec$Gene.Name
rownames(reads_empai) <- reads_spec$Gene.Name

reads_spec$accession <- NULL
reads_spec$Gene.Name <- NULL
reads_empai$accession <- NULL
reads_empai$Gene.Name <- NULL

reads_spec <- reads_spec %>% filter(!grepl("Reverse_", rownames(.)))
reads_empai <- reads_empai[rownames(reads_spec), , drop = FALSE]

feature_data <- full_data_df[, c("accession", "Gene.Name", "Description")]
feature_data <- feature_data[match(rownames(reads_spec), feature_data$Gene.Name), , drop = FALSE]
feature_data <- cbind(feature_data, extract_protein_fields(feature_data$Description))
rownames(feature_data) <- feature_data$Gene.Name

sample_index <- as.integer(sub("_.*$", "", sample_files))
metadata <- data.frame(
  Original.ID = sample_index,
  Sample.Names = sample_names,
  stringsAsFactors = FALSE
)

metadata <- tidyr::separate(metadata, Sample.Names, into = c("Particle", "Time", "Dx", "Sample.Order"), sep = "_", remove = FALSE)
metadata$Sample.Order <- ifelse(metadata$Sample.Order == "bis", metadata$Original.ID, metadata$Sample.Order)
metadata$Sample.Order <- as.numeric(metadata$Sample.Order)
metadata <- metadata[order(metadata$Sample.Order), ]
rownames(metadata) <- metadata$Sample.Names
metadata$short_setup <- paste0(metadata$Particle, "_", metadata$Time)

color_mapping <- c(
  "Au_4h" = "green",
  "PEG_4h" = "red",
  "Au_24h" = "darkgreen",
  "PEG_24h" = "firebrick4"
)
metadata$Fill <- "white"
metadata$Color <- unname(color_mapping[metadata$short_setup])

write_txt(reads_spec, file.path(paths$txt, "reads_spec.txt"))
write_txt(reads_empai, file.path(paths$txt, "reads_empai.txt"))
write_txt(feature_data, file.path(paths$txt, "feature_data.txt"))
write_txt(metadata, file.path(paths$txt, "metadata.txt"))

write_xlsx(reads_spec, file.path(paths$xlsx, "reads_spec.xlsx"))
write_xlsx(reads_empai, file.path(paths$xlsx, "reads_empai.xlsx"))
write_xlsx(feature_data, file.path(paths$xlsx, "feature_data.xlsx"))
write_xlsx(metadata, file.path(paths$xlsx, "metadata.xlsx"))
