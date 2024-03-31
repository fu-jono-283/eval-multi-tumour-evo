# Create the data for the table
data <- data.frame(
  Method = c("BulkSEQ", "scSEQ"),
  Advantages = c(
    "High throughput\nCost-effective for large sample sizes\nSuitable for studying average gene expression across cell populations",
    "Single-cell resolution\nProvides insights into cellular heterogeneity\nReveals rare cell populations"
  ),
  Disadvantages = c(
    "Limited resolution, averaging gene expression across cell populations\nMasking of rare cell populations\nCannot capture cellular heterogeneity",
    "High cost per cell\nIncreased technical variability\nSusceptible to amplification biases"
  )
)

# Split the advantages and disadvantages into separate rows
advantages_rows <- unlist(strsplit(data$Advantages, "\n"))
disadvantages_rows <- unlist(strsplit(data$Disadvantages, "\n"))

# Create a new data frame with separate rows for each advantage and disadvantage
new_data <- data.frame(
  Method = rep(data$Method, each = length(advantages_rows)),
  Category = rep(c("Advantages", "Disadvantages"), each = length(advantages_rows)),
  Point = c(advantages_rows, disadvantages_rows)
)

# Save the data as a CSV file
write.csv(new_data, "table_data.csv", row.names = FALSE)
