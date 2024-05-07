#CONSENSUS SEQUENCE CONCATENATION ANALYSIS ---------------------------------------------------------

library(ape)
library(TreeDist)
library(dplyr)


calculate_tree_distance <- function(tree1, tree2) {
    tree_distance <- TreeDistance(tree1, tree2)
    return(tree_distance)
}

csv_file_path <- "dataframes/testing_df/big_scores.csv"

combined_df <- read.csv(csv_file_path, stringsAsFactors = FALSE)

species_tree <- read.tree(Sys.getenv("species_tree_file"))
concat_tree <- read.tree(Sys.getenv("concat_besttree_file"))
population_size <- (Sys.getenv("population_size"))

tree_distance <- calculate_tree_distance(species_tree, concat_tree)

print(tree_distance)

combined_df <- combined_df %>%
    mutate(concat_seq_score = if_else(effective_population_size == population_size, tree_distance, concat_seq_score))

print(combined_df[["effective_population_size"]])
print(population_size)

write.csv(combined_df, csv_file_path, row.names = FALSE)



