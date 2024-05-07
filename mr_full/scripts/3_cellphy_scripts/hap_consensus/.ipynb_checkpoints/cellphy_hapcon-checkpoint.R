#HAPLOTYPE CONSENSUS SEQUENCE ANALYSIS ---------------------------------------------------------

library(ape)
library(TreeDist)
library(dplyr)

calculate_tree_distance <- function(tree1, tree2) {
    tree_distance <- TreeDistance(tree1, tree2)
    return(tree_distance)
}

iteration_identifier <- Sys.getenv("iteration_identifier")

csv_filename <- paste0("data/2_ccoal/", iteration_identifier, "/combined_data_200.csv")

species_tree <- read.tree(Sys.getenv("species_tree_file"))
hapcon_tree <- read.tree(Sys.getenv("hapcon_besttree_file"))
population_size <- (Sys.getenv("population_size"))
replicate <- Sys.getenv("replicate")

tree_distance <- calculate_tree_distance(species_tree, hapcon_tree)

combined_df <- read.csv(csv_filename, stringsAsFactors = FALSE)

combined_df$hapcon_seq_score <- ifelse(
  combined_df$effective_population_size == population_size & combined_df$Replicate == replicate,
  tree_distance,
  combined_df$hapcon_seq_score
)

write.csv(combined_df, csv_filename, row.names = FALSE)

cat("Score updated in the CSV file: ", csv_filename, "\n")
