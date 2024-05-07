source_file <- "dataframes/part1_template/combined_data_200.csv"

target_directory <- "data/2_ccoal"
subdirectories <- list.dirs(target_directory, full.names = TRUE, recursive = FALSE)

for (subdir in subdirectories) {
  target_file <- file.path(subdir, "combined_data_200.csv")
  if (!file.exists(target_file)) {
    file.copy(source_file, target_file, overwrite = TRUE)
  
    cat("Copied", source_file, "to", target_file, "\n")
  } else {
    cat("Skipped copying", source_file, "to", target_file, "because it already exists.\n")
  }
}