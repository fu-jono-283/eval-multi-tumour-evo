library(dplyr)

subdirs <- list.dirs("data/2_ccoal", full.names = TRUE, recursive = FALSE)

template_df <- read.csv("dataframes/part1_template/combined_data_200.csv")

compiled_hapcon_seq_scores <- rep(-1, nrow(template_df))
compiled_vcfbulk_seq_scores <- rep(-1, nrow(template_df))
compiled_single_seq_scores <- rep(-1, nrow(template_df))

compiled_hapcon_IUPAC_seq_scores <- rep(-1, nrow(template_df))

for (subdir in subdirs) {
  csv_file <- file.path(subdir, "combined_data_200.csv")
  
  if (file.exists(csv_file)) {
    df <- read.csv(csv_file)
    
    updated_hapcon_seq_scores <- df$hapcon_seq_score
    updated_vcfbulk_seq_scores <- df$vcfbulk_seq_score
    updated_single_seq_scores <- df$scSEQ_score
      
    updated_hapcon_IUPAC_seq_scores <- df$hapcon_IUPAC_seq_score  
    
    updated_hapcon_indices <- which(updated_hapcon_seq_scores >= 0)
    updated_vcfbulk_indices <- which(updated_vcfbulk_seq_scores >= 0)
    updated_single_indices <- which(updated_single_seq_scores >= 0)
      
    updated_hapcon_IUPAC_indices <- which(updated_hapcon_IUPAC_seq_scores >= 0)

    compiled_hapcon_seq_scores[updated_hapcon_indices] <- updated_hapcon_seq_scores[updated_hapcon_indices]
    compiled_vcfbulk_seq_scores[updated_vcfbulk_indices] <- updated_vcfbulk_seq_scores[updated_vcfbulk_indices]
    compiled_single_seq_scores[updated_single_indices] <- updated_single_seq_scores[updated_single_indices]
      
    compiled_hapcon_IUPAC_seq_scores[updated_hapcon_IUPAC_indices] <- updated_hapcon_IUPAC_seq_scores[updated_hapcon_IUPAC_indices]  
  }
}

template_df$vcfbulk_seq_score <- compiled_vcfbulk_seq_scores
template_df$scSEQ_score <- compiled_single_seq_scores
template_df$hapcon_seq_score <- compiled_hapcon_seq_scores

template_df$hapcon_IUPAC_seq_score <- compiled_hapcon_IUPAC_seq_scores

write.csv(template_df, "dataframes/testing_df/small_scores.csv", row.names = FALSE)
