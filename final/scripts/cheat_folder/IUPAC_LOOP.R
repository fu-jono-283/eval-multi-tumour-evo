base_dir <- "data/2_ccoal"
subdirectories <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

new_column_name <- "hapcon_IUPAC_seq_score"
new_column_data <- -1

for (subdir in subdirectories) {
  csv_file_path <- file.path(subdir, "combined_data_200.csv")

  if (file.exists(csv_file_path)) {
    combined_data <- read.csv(csv_file_path)
    combined_data[[new_column_name]] <- new_column_data  

    write.csv(combined_data, csv_file_path, row.names = FALSE)

    cat("Added", new_column_name, "to", csv_file_path, "\n")
  } else {
    cat("CSV file not found in", subdir, "\n")
  }
}
