create_dir <- function(path) {
  suppressWarnings(dir.create(path, recursive = TRUE, showWarnings = FALSE))
}

first_non_na <- function(values) {
  non_na <- values[!is.na(values) & values != ""]
  if (length(non_na) == 0) return(NA_character_)
  non_na[[1]]
}

as_numeric_safe <- function(x) {
  x <- gsub(';.*$', '', x)
  x <- gsub(';', '', x)
  x <- gsub('E', 'e', x)
  suppressWarnings(as.numeric(x))
}

build_threshold_tag <- function(threshold) {
  paste0('mean_gt_', threshold)
}
