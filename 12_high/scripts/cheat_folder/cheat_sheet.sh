find data/2_ccoal -type f -name "full_hap.0001_consensus.fasta" | wc -l

find data/2_ccoal/* -maxdepth 0 -type d '!' -exec test -e "{}/full_haplotypes_dir/full_hap.0001_consensus.fasta" ';' -print

find data/2_ccoal -mindepth 1 -maxdepth 1 -type d -exec cp {}/combined_data_120.csv {}/hapcon320edited.csv \;

find data/2_ccoal -type f -name "combined_data_120.csv" -exec sh -c 'cp "$0" "${0%combined_data_120.csv}hapcon320.csv"' {} \;

find data/2_ccoal -type f -name "hapcon320csv.csv" -exec sh -c 'mv "$0" "${0%hapcon320csv.csv}backupdf.csv"' {} \;

##########

#removing columns

library(dplyr)
df <- read.csv("dataframes/spare_combined_data_120.csv")

filtered_df <- df %>%
  filter(Replicate == "r1") %>%
  select(-Replicate, -vcfbulk_seq_score, -scSEQ_score, -hapcon_seq_score) %>%
  write.csv(file = "dataframes/spare_combined_data_120_filtered.csv", row.names = FALSE)

#########

#adding columns

csv_file <- "dataframes/spare_combined_data_120_filtered.csv"
df <- read.csv(csv_file)
df$concat_seq_score <- -1
df$ASTRAL_score <- -1
write.csv(df, file = csv_file, row.names = FALSE)
cat("Column 'concat_seq_score' added and CSV updated:", csv_file, "\n")

#########

#manually adding a score 

data <- read.csv("dataframes/combined_data_27.csv")

row_index <- 0
col_name <- "single_tree_scores"
new_value <- 0
data[row_index, col_name] <- new_value

write.csv(data, "dataframes/combined_data_27.csv", row.names = FALSE)
cat("Cell updated successfully.\n")

#########

#copying dataframe to every subfolder... 

source_file <- "dataframes/combined_data_120.csv"

target_directory <- "data/2_ccoal"
subdirectories <- list.dirs(target_directory, full.names = TRUE, recursive = FALSE)

for (subdir in subdirectories) {
  target_file <- file.path(subdir, "combined_data_120.csv")
  if (!file.exists(target_file)) {
    file.copy(source_file, target_file, overwrite = TRUE)
  
    cat("Copied", source_file, "to", target_file, "\n")
  } else {
    cat("Skipped copying", source_file, "to", target_file, "because it already exists.\n")
  }
}

#########

#merging the scores of all the subfolders' dataframe

library(dplyr)

subdirs <- list.dirs("data/2_ccoal", full.names = TRUE, recursive = FALSE)

template_df <- read.csv("dataframes/part1_template/combined_data_200.csv")

compiled_hapcon_seq_scores <- rep(-1, nrow(template_df))
compiled_vcfbulk_seq_scores <- rep(-1, nrow(template_df))
compiled_single_seq_scores <- rep(-1, nrow(template_df))

for (subdir in subdirs) {
  csv_file <- file.path(subdir, "combined_data_200.csv")
  
  if (file.exists(csv_file)) {
    df <- read.csv(csv_file)
    
    updated_hapcon_seq_scores <- df$hapcon_seq_score
    updated_vcfbulk_seq_scores <- df$vcfbulk_seq_score
    updated_single_seq_scores <- df$scSEQ_score
    
    updated_hapcon_indices <- which(updated_hapcon_seq_scores >= 0)
    updated_vcfbulk_indices <- which(updated_vcfbulk_seq_scores >= 0)
    updated_single_indices <- which(updated_single_seq_scores >= 0)

    compiled_hapcon_seq_scores[updated_hapcon_indices] <- updated_hapcon_seq_scores[updated_hapcon_indices]
    compiled_vcfbulk_seq_scores[updated_vcfbulk_indices] <- updated_vcfbulk_seq_scores[updated_vcfbulk_indices]
    compiled_single_seq_scores[updated_single_indices] <- updated_single_seq_scores[updated_single_indices]
  }
}

template_df$vcfbulk_seq_score <- compiled_vcfbulk_seq_scores
template_df$scSEQ_score <- compiled_single_seq_scores
template_df$hapcon_seq_score <- compiled_hapcon_seq_scores

write.csv(template_df, "dataframes/testing_df/small_scores.csv", row.names = FALSE)

#########

#this looks into the ccoal folder... and removes the log file???

main_folder <- "data/2_ccoal"

iteration_folders <- list.files(main_folder, full.names = TRUE)

for (iteration_folder in iteration_folders) {
 
  files_in_folder <- list.files(iteration_folder, full.names = TRUE)
  
  if (length(files_in_folder) == 1 && basename(files_in_folder) == "log") {
   
    unlink(iteration_folder, recursive = TRUE)
    cat("Removed folder:", iteration_folder, "\n")
  }
}

cat("Cleanup completed.\n")

#########

#this reads a csv, and removes the specified row 

data <- read.csv("dataframes/combined_data_160.csv")
data <- data[-1001, ]
write.csv(data, "dataframes/combined_data_160.csv", row.names = FALSE)
cat("Row 1001 has been removed from combined_data_160.csv\n")

#########

