list_excel_files <- function(path) {
  list.files(path = path, pattern = "\\.xlsx$", full.names = FALSE)
}

read_sample_table <- function(file_path) {
  as.data.frame(readxl::read_excel(file_path))
}

write_txt <- function(x, file_path, row_names = TRUE) {
  utils::write.table(x, file = file_path, row.names = row_names)
}

write_xlsx <- function(x, file_path) {
  openxlsx::write.xlsx(x, file = file_path, overwrite = TRUE)
}
