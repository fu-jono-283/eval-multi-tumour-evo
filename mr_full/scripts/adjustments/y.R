# I think this is updating scores manually... like outside of the main cellphy script

library(ape)
library(TreeDist)
library(dplyr)

calculate_tree_distance <- function(tree1, tree2) {
    tree_distance <- TreeDistance(tree1, tree2)
    return(tree_distance)
}

ccoal_parent_directory <- "data/2_ccoal"
ccoal_iteration_folders <- list.dirs(ccoal_parent_directory, full.names = FALSE, recursive = FALSE)

for (iteration_folder in ccoal_iteration_folders) {

    single_tree_path <- file.path(ccoal_parent_directory, iteration_folder, "single.raxml.bestTree")
    hapcon_tree_path <- file.path(ccoal_parent_directory, iteration_folder, "hapcon.raxml.bestTree")
    
    if (file.exists(hapcon_tree_path)) {
        hapcon_tree <- read.tree(hapcon_tree_path)
        
    population_size <- as.numeric(gsub(".*_pop(\\d+)_sites\\d+r\\d+", "\\1", iteration_folder))
    replicate <- gsub(".*?(r\\d+)$", "\\1", iteration_folder)

    directory_path <- file.path("data/1_simphy", iteration_folder)
        
    species_tree_file <- list.files(directory_path, pattern = "s_tree.trees", full.names = TRUE)    
    gene_tree_file <- list.files(directory_path, pattern = "^g_trees\\d+\\.trees$", full.names = TRUE)

    species_tree <- read.tree(species_tree_file)    
    gene_tree <- read.tree(gene_tree_file)

    tree_distance <- calculate_tree_distance(species_tree, hapcon_tree)
    
    csv_filename <- file.path(ccoal_parent_directory, iteration_folder, "combined_data_200.csv")

    combined_df <- read.csv(csv_filename, stringsAsFactors = FALSE)

combined_df$hapcon_seq_score <- ifelse(
  combined_df$effective_population_size == population_size & combined_df$Replicate == replicate,
  tree_distance,
  combined_df$hapcon_seq_score
)
print(population_size)
print(replicate)

    write.csv(combined_df, csv_filename, row.names = FALSE)

    cat("Score updated in the CSV file for ccoal iteration: ", csv_filename, "\n")
}
}