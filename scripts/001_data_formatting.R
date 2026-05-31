suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(openxlsx)
  library(readxl)
})

source("src/utils.R")
source("src/io.R")
source("src/preprocessing.R")

# Step 1: load all raw input files and standardize data types/content.
paths <- list(
  raw = "Inputs/001_Raw_data",
  quality = "Inputs/002_Quality_annotated",
  txt = "Inputs/txt",
  xlsx = "Inputs/xlsx"
)

create_dir(paths$quality)
create_dir(paths$txt)
create_dir(paths$xlsx)

sample_names <- list_excel_files(paths$raw)

dim_results_df <- data.frame(Sample.Name = character(), N.genes = integer(), stringsAsFactors = FALSE)
dupe_gene_name_list <- list()
dupe_row_list <- list()

for (sample in sample_names) {
  message("Processing: ", sample)
  tab <- read_sample_table(file.path(paths$raw, sample))
  tab <- clean_sample_table(tab)

  protein_fields <- extract_protein_fields(tab$description)
  tab <- cbind(tab, protein_fields)

  write_xlsx(tab, file.path(paths$quality, sample))

  dim_results_df <- rbind(dim_results_df, data.frame(Sample.Name = sample, N.genes = nrow(tab)))

  dupe_gene_name <- tab[tab$Gene.Name.issues == "Duplicates", ]
  dupe_gene_name_list[[sample]] <- dupe_gene_name

  dupe_row <- tab[tab$duplicated_row == "Duplicates", ]
  dupe_row_list[[sample]] <- dupe_row
}

dupe_gene_name_list_filtered <- Filter(function(x) nrow(x) > 0, dupe_gene_name_list)
dupe_row_list_filtered <- Filter(function(x) nrow(x) > 0, dupe_row_list)

write_txt(dim_results_df, file.path(paths$txt, "Num_genes_per_condition.txt"), row_names = FALSE)
write_xlsx(dim_results_df, file.path(paths$xlsx, "Num_genes_per_condition.xlsx"))
write_xlsx(dupe_gene_name_list_filtered, file.path(paths$xlsx, "dupe_gene_name_condition.xlsx"))
write_xlsx(dupe_row_list_filtered, file.path(paths$xlsx, "dupe_rows_condition.xlsx"))
