library(ape)
library(TreeDist)
library(dplyr)

source_directory <- "data/2_ccoal"
destination_directory <- "data/9_true_single"

df <- read.csv("dataframes/testing_df/big_scores-Copy1.csv")
df$CCOAL_ASTRAL_score <- -1
write.csv(df, file = "dataframes/testing_df/big_scores-Copy1.csv", row.names = FALSE)

process_tree_files <- function(input_file_path, output_file_path) {
    tree_content <- readLines(input_file_path, warn = FALSE)
    tree_content <- sub("^\\(", "", tree_content)
    tree_content <- gsub(":\\d+\\.\\d+,outgcell:\\d+\\.\\d+\\)", "", tree_content)
    writeLines(tree_content, output_file_path)
  }

replicate_directories <- list.files(source_directory, full.names = TRUE)

for (replicate_dir in replicate_directories) {
  dir_name <- basename(replicate_dir)
  pop_sites <- gsub("r[0-9]+$", "", dir_name)
  
  dest_subdir <- file.path(destination_directory, pop_sites)
  dir.create(dest_subdir, recursive = TRUE, showWarnings = FALSE)

  input_file_path <- file.path(replicate_dir, "trees_dir")
  output_file_path <- file.path(replicate_dir, "trees_dir")
    
  tree_proc <- list.files(input_file_path, pattern = "trees.0001", full.names = TRUE)

  for (input_file_path in tree_proc) {
    output_file_name <- paste0(tools::file_path_sans_ext(basename(input_file_path)), "_outgcell_removed")
    output_file_path <- file.path(output_file_path, output_file_name)
    process_tree_files(input_file_path, output_file_path)
  
  }

tree_files <- list.files(file.path(replicate_dir, "trees_dir"), pattern = "trees_outgcell_removed$", full.names = TRUE)

    
    print(tree_files)

  for (tree_file in tree_files) {
    replicate_label <- gsub(".*r([0-9]+)$", "r\\1", replicate_dir)
    dest_file <- file.path(dest_subdir, paste0(replicate_label, "_", basename(tree_file)))
    file.copy(tree_file, dest_file)
      
  }  
}

python_script <- "scripts/9_true_single_scripts/merge_trees.py"
system2("python", args = python_script)



python_script_2 <- "scripts/9_true_single_scripts/mapping_format.py"
system2("python", args = python_script_2)




iteration_folders <- list.dirs("data/9_true_single", full.names = TRUE, recursive = FALSE)

for (iteration_folder in iteration_folders) {
 if (iteration_folder != "data/9_true_single" && !grepl("/\\.ipynb_checkpoints$", iteration_folder)) {
pop_sites_match <- regmatches(iteration_folder, gregexpr("_pop[0-9]+_sites[0-9]+", iteration_folder))
 pop_sites <- pop_sites_match[[1]]
      
  simphy_replicate_1_folder <- file.path("data/1_simphy", paste0(pop_sites, "r1"))
  species_tree_file <- file.path(simphy_replicate_1_folder, "s_tree.trees")
  species_tree <- file.path(iteration_folder, "s_tree.trees")
    
  cat("Copying species tree from:", species_tree_file, "to:", species_tree, "\n") 
  file.copy(species_tree_file, species_tree)   

  mapping_file <- file.path(iteration_folder, "mapping.txt")
  merged_trees_file <- file.path(iteration_folder, "merged_trees.tre")
  astral_tree_file <- file.path(iteration_folder, "astral_tree")

  astral_command <- paste("astral -a", mapping_file, "-i", merged_trees_file, "-o", astral_tree_file)
  system(astral_command)   
}}

csv_filename <- "dataframes/testing_df/big_scores-Copy1.csv"
combined_df <- read.csv(csv_filename, stringsAsFactors = FALSE)

for (iteration_folder in iteration_folders) {
  astral_tree_file <- file.path(iteration_folder, "astral_tree")
  species_tree_file <- file.path(iteration_folder, "s_tree.trees")
  
  if (file.exists(astral_tree_file) && file.size(astral_tree_file) > 0 &&
      file.exists(species_tree_file) && file.size(species_tree_file) > 0) {
    astral_tree <- read.tree(astral_tree_file)
    species_tree <- read.tree(species_tree_file)
      
print(astral_tree)
print(species_tree)
      
    calculate_tree_distance <- function(tree1, tree2) {
      tree_distance <- TreeDistance(tree1, tree2)
        print(tree_distance)
      return(tree_distance)
    }
    
    tree_distance <- calculate_tree_distance(astral_tree, species_tree)      
    population_size <- as.numeric(gsub(".*_pop([0-9]+)_sites([0-9]+).*", "\\1", iteration_folder))

  combined_df <- combined_df %>%
    mutate(
      CCOAL_ASTRAL_score = ifelse(
        effective_population_size == population_size,
        tree_distance,
        CCOAL_ASTRAL_score
      )
    )
      
   write.csv(combined_df, file = "dataframes/testing_df/big_scores.csv", row.names = FALSE)
  cat("CSV file updated.\n")
}
}
