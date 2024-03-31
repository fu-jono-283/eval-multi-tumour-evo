library(ape)
library(TreeDist)
library(dplyr)

source_directory <- "data/2_ccoal"
destination_directory <- "data/6_super"

replicate_directories <- list.files(source_directory, full.names = TRUE)

for (replicate_dir in replicate_directories) {
  dir_name <- basename(replicate_dir)
  pop_sites <- gsub("r[0-9]+$", "", dir_name)
  
  dest_subdir <- file.path(destination_directory, pop_sites)
  dir.create(dest_subdir, recursive = TRUE, showWarnings = FALSE)
    
  tree_files <- list.files(replicate_dir, pattern = "hapcon.raxml.bestTree", full.names = TRUE)
  
  for (tree_file in tree_files) {
    replicate_label <- gsub(".*r([0-9]+)$", "r\\1", replicate_dir)
    dest_file <- file.path(dest_subdir, paste0(replicate_label, "_", basename(tree_file)))
    file.copy(tree_file, dest_file)
  }  
}

python_script <- "scripts/6_supertree_scripts/merge_trees.py"
system2("python", args = python_script)

iteration_folders <- list.dirs("data/6_super", full.names = TRUE, recursive = FALSE)

for (iteration_folder in iteration_folders) {
 if (iteration_folder != "data/6_super" && !grepl("/\\.ipynb_checkpoints$", iteration_folder)) {
pop_sites_match <- regmatches(iteration_folder, gregexpr("_pop[0-9]+_sites[0-9]+", iteration_folder))
 pop_sites <- pop_sites_match[[1]]
      
  simphy_replicate_1_folder <- file.path("data/1_simphy", paste0(pop_sites, "r1"))
  species_tree_file <- file.path(simphy_replicate_1_folder, "s_tree.trees")
  species_tree <- file.path(iteration_folder, "s_tree.trees")
    
  cat("Copying species tree from:", species_tree_file, "to:", species_tree, "\n") 
  file.copy(species_tree_file, species_tree)   

  merged_trees_file <- file.path(iteration_folder, "merged_trees.tre")
  supertree_file <- file.path(iteration_folder, "supertree")
     
print(merged_trees_file)
print(supertree_file)

  fastRFS_command <- paste("FastRFS -i", merged_trees_file, "-o", supertree_file)
  system(fastRFS_command)
}}

csv_filename <- "dataframes/testing_df/big_scores.csv"
combined_df <- read.csv(csv_filename, stringsAsFactors = FALSE)

for (iteration_folder in iteration_folders) {
  supertree_file <- file.path(iteration_folder, "supertree.single")
  species_tree_file <- file.path(iteration_folder, "s_tree.trees")
  
  if (file.exists(supertree_file) && file.size(supertree_file) > 0 &&
      file.exists(species_tree_file) && file.size(species_tree_file) > 0) {
    supertree <- read.tree(supertree_file)
    species_tree <- read.tree(species_tree_file)
      
print(supertree)
print(species_tree)
      
    calculate_tree_distance <- function(tree1, tree2) {
      tree_distance <- TreeDistance(tree1, tree2)
        print(tree_distance)
      return(tree_distance)
    }
    
    tree_distance <- calculate_tree_distance(supertree, species_tree)      
    population_size <- as.numeric(gsub(".*_pop([0-9]+)_sites([0-9]+).*", "\\1", iteration_folder))

  combined_df <- combined_df %>%
    mutate(
      hap_FastRFS_score = ifelse(
        effective_population_size == population_size,
        tree_distance,
        hap_FastRFS_score
      )
    )
      
   write.csv(combined_df, file = csv_filename, row.names = FALSE)
  cat("CSV file updated.\n")
}
}





