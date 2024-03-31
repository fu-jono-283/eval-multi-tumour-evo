library(ape)
library(TreeDist)
library(dplyr)

source_directory <- "data/2_ccoal"
destination_directory <- "data/4_astral"

df <- read.csv("dataframes/part1_template/combined_data.csv")

filtered_df <- df %>%
  filter(Replicate == "r1") %>%
  select(-Replicate, -scSEQ_score, -hapconsensus_score, -vcfpseudobulk_score) %>%
  write.csv(file = "dataframes/part2_template/combined_data_BIG.csv", row.names = FALSE)

# Constructing second dataframe -----------------------------------------------------------------------------------

df2 <- read.csv("dataframes/part2_template/combined_data_BIG.csv")
df2$ASTRAL_score <- -1
df2$concat_hap_score <- -1 
df2$FastRFS_hap_score <- -1
df2$concat_vcf_score <- -1 
df2$FastRFS_vcf_score <- -1
write.csv(df2, file = "dataframes/part2_template/combined_data_BIG.csv", row.names = FALSE)

replicate_directories <- list.files(source_directory, full.names = TRUE)

for (replicate_dir in replicate_directories) {
  dir_name <- basename(replicate_dir)
  pop_sites <- gsub("r[0-9]+$", "", dir_name)
   
  dest_subdir <- file.path(destination_directory, pop_sites)
  dir.create(dest_subdir, recursive = TRUE, showWarnings = FALSE)
    
  tree_files <- list.files(replicate_dir, pattern = "single.raxml.bestTree", full.names = TRUE)
  
    
  for (tree_file in tree_files) {
    replicate_label <- gsub(".*r([0-9]+)$", "r\\1", replicate_dir)
    dest_file <- file.path(dest_subdir, paste0(replicate_label, "_", basename(tree_file)))
    file.copy(tree_file, dest_file)
  }  
}

# ASTRAL--------------------------------------------------------------------------------------------------

python_script <- "scripts/4_astral_scripts/astral_merge_trees.py"
system2("python", args = python_script)

python_script_2 <- "scripts/4_astral_scripts/astral_mapping_format.py"
system2("python", args = python_script_2)

iteration_folders <- list.dirs("data/4_astral", full.names = TRUE, recursive = FALSE)

for (iteration_folder in iteration_folders) {
 if (iteration_folder != "data/4_astral" && !grepl("/\\.ipynb_checkpoints$", iteration_folder)) {
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
  print(astral_command)
  system(astral_command)   
}}

csv_filename <- "dataframes/part2_template/combined_data_BIG.csv"
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
      ASTRAL_score = ifelse(
        effective_population_size == population_size,
        tree_distance,
        ASTRAL_score
      )
    )
      
   write.csv(combined_df, file = "dataframes/testing_df/big_scores.csv", row.names = FALSE)
  cat("CSV file updated.\n")
}
}

# CONCATENATED HAPLOTYPE CONSENSUS SEQUENCES--------------------------------------------------------------------------------------------------

initial_dir <- "data/2_ccoal"
main_dir <- "data/5_concat_hap/"

subdirs <- list.dirs(initial_dir, full.names = TRUE, recursive = TRUE) 
subdirs <- subdirs[!grepl("_ipynb_checkpoints", subdirs)]

for (subdir in subdirs) {
  consensus_file <- file.path(subdir, "full_haplotypes_dir", "full_hap.0001_consensus.fasta")
  if (file.exists(consensus_file)) {
    replicate_label <- basename(subdir)
    target_dir <- file.path(main_dir, replicate_label)
    if (!dir.exists(target_dir)) {
      dir.create(target_dir, recursive = TRUE)
    }
    file.copy(consensus_file, target_dir)
  }
}

main_subdirs <- list.dirs(main_dir, full.names = TRUE, recursive = FALSE)
subdir_identifiers <- list()

for (subdir in main_subdirs) {
  identifier <- gsub(".*(_pop\\d+_sites\\d+r)\\d+.*", "\\1", subdir)
    if (!(identifier %in% names(subdir_identifiers))) {
    subdir_identifiers[[identifier]] <- list()
  }
    subdir_identifiers[[identifier]] <- c(subdir_identifiers[[identifier]], subdir)
}

for (identifier in names(subdir_identifiers)) {
  subdir_w_identifiers <- subdir_identifiers[[identifier]]
     merged_folder <- file.path(main_dir, gsub("\\W", "", identifier))
    if (!file.exists(merged_folder)) {
    dir.create(merged_folder)
  }
  
  replicate_counter <- 1
  
  for (folder in subdir_w_identifiers) {
    files <- list.files(folder)
    for (file in files) {
      if (grepl("^full_hap\\.\\d+_consensus\\.fasta$", file)) {
        current_file <- file.path(folder, file)
        replicate_label <- sprintf("r%d", replicate_counter)
        new_file <- file.path(merged_folder, paste0("full_hap.0001_consensus_", replicate_label, ".fasta"))
        file.rename(current_file, new_file)
        cat("Moved", current_file, "to", new_file, "\n")
        replicate_counter <- replicate_counter + 1
      }
    }    
    if (file.info(folder)$isdir) {
      unlink(folder, recursive = TRUE)
      cat("Removed folder:", folder, "\n")
    }
  }
    
source_species_tree <- file.path("data/1_simphy", paste0(identifier, "1"), "s_tree.trees")
print(paste("Debug - Source Species Tree for Characteristic", identifier, ":", source_species_tree))

target_species_tree <- file.path(merged_folder, "s_tree.trees")

  if (file.exists(source_species_tree)) {
    file.copy(source_species_tree, target_species_tree)
    cat("Species tree copied to:", target_species_tree, "\n")
  } else {
    cat("Error: Species tree not found for identifier", identifier, "\n")
  }
}

concat_dir <- "data/5_concat_hap/"
subdirs <- list.dirs(concat_dir, full.names = TRUE, recursive = FALSE)

for (subdir in subdirs) {
  if (grepl("r$", subdir)) {
    new_subdir <- gsub("r$", "", subdir)  
    if (subdir != new_subdir) {
      file.rename(subdir, new_subdir)
      cat("Renamed folder:", subdir, "to", new_subdir, "\n")
    }
  }
}

main_folder <- "data/5_concat_hap"

subdirectories <- list.dirs(main_folder, full.names = TRUE, recursive = FALSE)

sequences <- list()

for (folder_path in subdirectories) {
  sequences <- list()
  
  fasta_files <- list.files(folder_path, pattern = ".fasta$", full.names = TRUE)
  
  for (file in fasta_files) {
    lines <- readLines(file)
    
    current_header <- NULL
    current_sequence <- NULL
    
    for (line in lines) {
      if (startsWith(line, ">")) {
        if (!is.null(current_header)) {
          sequences[[current_header]] <- paste(sequences[[current_header]], current_sequence, sep = "")
        }
        current_header <- line
        current_sequence <- ""
      } else {
        current_sequence <- paste(current_sequence, line, sep = "")
      }
    }
    
    sequences[[current_header]] <- paste(sequences[[current_header]], current_sequence, sep = "")
  }
  
  output_file <- file.path(folder_path, "concat.fasta")
  
  writeLines(paste(names(sequences), unlist(sequences), sep = "\n"), con = output_file)
}

#CellPhy for HCS 

check_condition <- function(bestTree_files) {
  return(length(bestTree_files) == 8) # adjust 
}

script_folder <- "scripts/hap_concat"
script_pattern <- "^pre.*\\.R$"
script_files <- list.files(path = script_folder, pattern = script_pattern, full.names = TRUE)

while (TRUE) {
  for (script_file in script_files) {
    cat("Executing script:", script_file, "\n")
    system(paste("Rscript", shQuote(script_file)))
  }

  files <- list.files("data/2_ccoal", pattern = "^hapcon\\.raxml\\.bestTree$", full.names = TRUE, recursive = TRUE)

  cat("List of best tree files:", files, "\n")

  condition_met <- check_condition(files)
  if (condition_met) {
    cat("All trees have been produced. Exiting the loop.\n")
    break
  }
  Sys.sleep(300)
}



