data <- read.csv("dataframes/testing_df/big_scores.csv")

row_index <- 145
col_name <- "concat_seq_score"
new_value <- 0
data[row_index, col_name] <- new_value

write.csv(data, "dataframes/testing_df/big_scores.csv", row.names = FALSE)
cat("Cell updated successfully.\n")