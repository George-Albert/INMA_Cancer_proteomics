clean_sample_table <- function(tab) {
  semicolon_rows <- grep("(.+);(.+)", tab$NSAF)
  semicolon_rows_empai <- grep("(.+);(.+)", tab$EMPAI)

  if (length(semicolon_rows) > 0) {
    tab[semicolon_rows, "NSAF"] <- stringr::str_replace_all(tab[semicolon_rows, "NSAF"], ";(.+)", "")
  }
  if (length(semicolon_rows_empai) > 0) {
    tab[semicolon_rows_empai, "EMPAI"] <- stringr::str_replace_all(tab[semicolon_rows_empai, "EMPAI"], ";(.+)", "")
  }

  for (col in seq_along(tab)) {
    if (is.character(tab[, col])) {
      tab[, col] <- stringr::str_replace_all(tab[, col], ";", "")
      tab[, col] <- stringr::str_replace_all(tab[, col], "E", "e")
    }
  }

  tab$seq.coverage <- as_numeric_safe(tab$seq.coverage)
  tab$NSAF <- as_numeric_safe(tab$NSAF)
  tab$EMPAI <- as_numeric_safe(tab$EMPAI)

  index_gene_name <- which(is.na(tab$Gene.Name))
  if (length(index_gene_name) > 0) {
    tab[index_gene_name, "Gene.Name"] <- tab[index_gene_name, "accession"]
  }

  tab$Gene.Name.issues <- "No"
  if (length(index_gene_name) > 0) {
    tab[index_gene_name, "Gene.Name.issues"] <- "contaminant"
  }

  duplicate_index <- which(duplicated(tab$Gene.Name))
  if (length(duplicate_index) > 0) {
    tab[duplicate_index, "Gene.Name.issues"] <- "Duplicates"
  }

  tab <- tab |>
    dplyr::mutate(
      duplicated_row = duplicated(tab[, c("seq.coverage", "seq.count", "spec.count", "NSAF", "EMPAI")]) |
        duplicated(tab[, c("seq.coverage", "seq.count", "spec.count", "NSAF", "EMPAI")], fromLast = TRUE)
    )

  tab$duplicated_row <- ifelse(tab$duplicated_row, "Duplicates", "No")

  remove_index <- grep(pattern = "contaminant", tab$accession)
  if (length(remove_index) > 0) {
    tab <- tab[-remove_index, ]
  }

  tab
}

extract_protein_fields <- function(description_vector) {
  data.frame(
    Protein.Name = gsub("^(.*?)\\sOS=.*$", "\\1", description_vector),
    Organism.Source = stringr::str_match(description_vector, "OS=(.*?)\\sOX=")[, 2],
    Taxonomy = stringr::str_match(description_vector, "OX=(.*?)\\s(GN=|Pe=)")[, 2],
    Gene.Name.description = stringr::str_match(description_vector, "GN=(.*?)\\sPe=")[, 2],
    Pe = stringr::str_match(description_vector, "Pe=(.*?)\\sSV=")[, 2],
    SV = stringr::str_match(description_vector, "SV=(\\d+)")[, 2]
  )
}

build_sample_variable_name <- function(file_name) {
  base <- tools::file_path_sans_ext(basename(file_name))
  sub("^([0-9]+)_(.*)$", "\\2_\\1", base)
}
