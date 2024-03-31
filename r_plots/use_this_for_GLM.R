# Define parameter labels
parameter_labels <- data.frame(
  Parameter = c("(Intercept)", "no_of_taxa", "exponential_growth_rate",
                "individuals_per_taxa", "effective_population_size",
                "species_birth_rate", "tree_wide_sub_rate", "no_of_sites",
                "ado", "amplification_error", "sequencing_errors"),
  Label = c("Intercept","Number of Tumours", "Exponential Growth Rate",
            "Cells per Tumour", "Effective Population Size",
            "Metastatic Rate", "Mutation Rate", "Number of Sites",
            "Allelic Dropout", "Amplification Error", "Sequencing Errors")
)

# Cleaning data 
data <- small_scores_low_SMALL
data$individuals_per_taxa <- as.numeric(data$individuals_per_taxa)
data <- na.omit(data)

# Function to perform GLM analysis
perform_glm <- function(data, score_column) {
  data[[score_column]] <- 1 - data[[score_column]]
  single_data <- data[data[[score_column]] != 2,] 
  
  glm_result <- glm(paste(score_column, "~",
                          "no_of_taxa + exponential_growth_rate + individuals_per_taxa +",
                          "effective_population_size + species_birth_rate + tree_wide_sub_rate +",
                          "no_of_sites + ado + amplification_error + sequencing_errors"), 
                    data = single_data)
  
  return(summary(glm_result))
}

# Perform GLM analysis for each score column
scSEQ_summary <- perform_glm(data, "scSEQ_score")
hapcon_summary <- perform_glm(data, "hapcon_seq_score")
vcfbulk_summary <- perform_glm(data, "vcfbulk_seq_score")

scSEQ_summary

# Function to add significance stars to coefficients
add_significance_stars <- function(summary_df) {
  stars <- ifelse(summary_df[, 4] < 0.001, "***",
                  ifelse(summary_df[, 4] < 0.01, "**",
                         ifelse(summary_df[, 4] < 0.05, "*",
                                "")))
  summary_df <- cbind(summary_df, stars)
  return(summary_df)
}

# Perform GLM analysis for each score column
perform_glm <- function(data, score_column) {
  data[[score_column]] <- 1 - data[[score_column]]
  single_data <- data[data[[score_column]] != 2,] 
  
  glm_result <- glm(paste(score_column, "~",
                          "no_of_taxa + exponential_growth_rate + individuals_per_taxa +",
                          "effective_population_size + species_birth_rate + tree_wide_sub_rate +",
                          "no_of_sites + ado + amplification_error + sequencing_errors"), 
                    data = single_data)
  
  # Extract coefficients from the summary object
  coefficients <- summary(glm_result)$coefficients
  
  # Add significance stars
  coefficients_with_stars <- add_significance_stars(coefficients)
  
  return(coefficients_with_stars)
}

# Add significance stars to each set of coefficients
scSEQ_summary_with_stars <- perform_glm(data, "scSEQ_score")
hapcon_summary_with_stars <- perform_glm(data, "hapcon_seq_score")
vcfbulk_summary_with_stars <- perform_glm(data, "vcfbulk_seq_score")

# Convert to data frames
scSEQ_df <- as.data.frame(scSEQ_summary_with_stars)
hapcon_df <- as.data.frame(hapcon_summary_with_stars)
vcfbulk_df <- as.data.frame(vcfbulk_summary_with_stars)

# Combine data frames with labels and stars
scSEQ_df <- cbind(Method = "scSEQ", scSEQ_df)
hapcon_df <- cbind(Method = "Hapcon", hapcon_df)
vcfbulk_df <- cbind(Method = "VCFBulk", vcfbulk_df)

# Combine data frames vertically
combined_df <- rbind(scSEQ_df, hapcon_df, vcfbulk_df)

# Create a column for parameter labels
combined_df$Parameter <- parameter_labels$Label

# Move the "Parameter" column to the second position
combined_df <- combined_df[, c(1, ncol(combined_df), 2:(ncol(combined_df)-1))]

# Reset row names to NULL
row.names(combined_df) <- NULL

# Print the dataframe
print(combined_df)
