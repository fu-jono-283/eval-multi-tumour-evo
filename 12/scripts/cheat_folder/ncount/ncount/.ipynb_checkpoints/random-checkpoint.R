# Read the CSV file back into R
ncount_data <- read.csv("dataframes/random_experiments/ncountdata.csv")

# Calculate total sites and add it as a new column
ncount_data$total_sites <- ncount_data$no_of_taxa * ncount_data$individuals_per_taxa * ncount_data$no_of_sites

# Calculate Proportion_N
ncount_data$Proportion_N <- ncount_data$N_count / ncount_data$total_sites

# Write the updated data to the CSV file
write.csv(ncount_data, file = "dataframes/random_experiments/ncountdata1.csv", row.names = FALSE)
