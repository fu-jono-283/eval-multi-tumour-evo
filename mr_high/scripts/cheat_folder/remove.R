data <- read.csv("dataframes/testing_df/big_scores.csv")
data <- data[-201, ]
write.csv(data, "dataframes/testing_df/big_scores.csv", row.names = FALSE)
cat("Row 1001 has been removed from combined_data_160.csv\n")