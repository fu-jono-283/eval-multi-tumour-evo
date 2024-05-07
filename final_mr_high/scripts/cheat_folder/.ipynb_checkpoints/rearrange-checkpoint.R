csv_file <- "dataframes/part1_template/combined_data_200.csv"
df <- read.csv(csv_file)

move_to_back <- df[c("hapcon_seq_score", "hapcon_IUPAC_seq_score")]

df <- df[, !(names(df) %in% c("hapcon_seq_score", "hapcon_IUPAC_seq_score"))]

df <- cbind(df, move_to_back)

write.csv(df, file = csv_file, row.names = FALSE)

cat("Columns moved to the back of the CSV:", csv_file, "\n")
