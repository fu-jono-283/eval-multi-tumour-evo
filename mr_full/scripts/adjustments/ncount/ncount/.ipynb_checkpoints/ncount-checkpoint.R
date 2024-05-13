library(Biostrings)
library(dplyr)

count_N_bases <- function(consensus_file) {
  records <- readDNAStringSet(consensus_file, format = "fasta")
  total_N_bases <- sum(vcountPattern("N" , records))
  return(total_N_bases)
}

single_ccoal_output_folder <- "data/2_ccoal"

cutoff_values <- seq(0, 1, by = 0.05)
print(cutoff_values)

combined_data <- data.frame(Iteration = numeric(),
                            CutoffValue = numeric(),
                            N_count = numeric(),
                            stringsAsFactors = FALSE)

iteration_folders <- list.files(path = single_ccoal_output_folder, full.names = TRUE)

for (iteration_folder in iteration_folders) {
  iteration_number <- as.numeric(gsub("^.*_pop(\\d+).*", "\\1", basename(iteration_folder)))
  replicate_label <- gsub("^.*r(\\d+)$", "r\\1", basename(iteration_folder))

  cutoff_counts <- list()
  
  for (cutoff in cutoff_values) {
    formatted_cutoff <- sprintf("%.2f", cutoff)
    full_hap_dir <- file.path(iteration_folder, "full_haplotypes_dir")
    consensus_file <- file.path(full_hap_dir, paste0("full_hap_consensus_",  formatted_cutoff, ".fasta"))
    cat("Checking consensus file:", consensus_file, "\n")

    if (file.exists(consensus_file)) {
      total_N_bases <- count_N_bases(consensus_file)
      cutoff_counts[[as.character(cutoff)]] <- total_N_bases
    }
  }
  

  if (length(cutoff_counts) > 0) {
    row <- data.frame(
      Iteration = iteration_number,
      CutoffValue = names(cutoff_counts),
      N_count = unlist(cutoff_counts),
      Replicate = replicate_label
    )
    combined_data <- rbind(combined_data, row)
  write.csv(combined_data, file = "dataframes/combined_data", row.names = FALSE)
  }
}

existing_data <- read.csv("dataframes/part1_template/combined_data_200.csv")
existing_data$CutoffValue <- NA
existing_data$N_count <- NA
existing_data$Ne <- existing_data$effective_population_size
existing_data$replabel <- existing_data$Replicate

mega_data <- merge(
  combined_data,
  existing_data,
  by.x = c("Iteration", "Replicate"),  
  by.y = c("Ne", "replabel"),         
  all.x = TRUE
)


mega_data$CutoffValue <- ifelse(is.na(mega_data$CutoffValue.x), mega_data$CutoffValue.y, mega_data$CutoffValue.x)
mega_data$N_count <- ifelse(is.na(mega_data$N_count.x), mega_data$N_count.y, mega_data$N_count.x)

mega_data <- mega_data[, !names(mega_data) %in% c("CutoffValue.x", "CutoffValue.y", "N_count.x", "N_count.y")]
mega_data <- mega_data[, c("CutoffValue", "N_count", setdiff(names(mega_data), c("Iteration", "CutoffValue", "N_count")))]

mega_data$Iteration <- NULL
mega_data$Proportion_N <- mega_data$N_count / mega_data$no_of_sites
mega_data <- subset(mega_data, Replicate == "r1")
mega_data <- mega_data %>%
  select(Proportion_N, everything())

write.csv(mega_data, file = "dataframes/random_experiments/ncountdata.csv", row.names = FALSE)



