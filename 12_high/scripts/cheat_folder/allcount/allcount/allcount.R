library(Biostrings)

count_ALL_bases <- function(consensus_file) {
  records <- readDNAStringSet(consensus_file, format = "fasta")
  total_ALL_bases <- sum(vcountPattern("N", records) +
                         vcountPattern("M", records) +
                         vcountPattern("K", records) +
                         vcountPattern("S", records) +
                         vcountPattern("W", records) +
                         vcountPattern("R", records) +
                         vcountPattern("Y", records))
  return(total_ALL_bases)
}


single_ccoal_output_folder <- "data/2_ccoal"

cutoff_values <- seq(0, 1, by = 0.05)
print(cutoff_values)

combined_data <- data.frame(Iteration = numeric(),
                            CutoffValue = numeric(),
                            ALL_count = numeric(),
                            stringsAsFactors = FALSE)

iteration_folders <- list.files(path = single_ccoal_output_folder, full.names = TRUE)

for (iteration_folder in iteration_folders) {
  iteration_number <- as.numeric(gsub("^.*_pop(\\d+).*", "\\1", basename(iteration_folder)))
  
  for (cutoff in cutoff_values) {
    formatted_cutoff <- sprintf("%.2f", cutoff)
    full_hap_dir <- file.path(iteration_folder, "full_haplotypes_dir")
    consensus_file <- file.path(full_hap_dir, paste0("full_hap_consensus_",  formatted_cutoff, ".fasta"))
    cat("Checking consensus file:", consensus_file, "\n")

    if (file.exists(consensus_file)) {
      total_ALL_bases <- count_ALL_bases(consensus_file)
      row <- data.frame(Iteration = iteration_number, CutoffValue = cutoff, ALL_count = total_ALL_bases)
      combined_data <- rbind(combined_data, row)
    }
  }
}

existing_data <- read.csv("dataframes/combined_data_120-Copy1.csv")
existing_data$CutoffValue <- NA
existing_data$ALL_count <- NA

mega_data <- merge(combined_data, existing_data, by.x = "Iteration", by.y = "effective_population_size", all.x = TRUE)
mega_data$CutoffValue <- ifelse(is.na(mega_data$CutoffValue.x), mega_data$CutoffValue.y, mega_data$CutoffValue.x)
mega_data$ALL_count <- ifelse(is.na(mega_data$ALL_count.x), mega_data$ALL_count.y, mega_data$ALL_count.x)
mega_data <- mega_data[, !names(mega_data) %in% c("CutoffValue.x", "CutoffValue.y", "ALL_count.x", "ALL_count.y")]
mega_data <- mega_data[, c("CutoffValue", "ALL_count", setdiff(names(mega_data), c("Iteration", "CutoffValue", "ALL_count")))]
mega_data$Iteration <- NULL
mega_data$Proportion_ALL <- mega_data$ALL_count / mega_data$no_of_sites
mega_data <- data.frame(Proportion_ALL = mega_data$Proportion_ALL, mega_data[, setdiff(names(mega_data), "Proportion_ALL")])

write.csv(mega_data, file = "dataframes/random_experiments/allcountdata.csv", row.names = FALSE)
