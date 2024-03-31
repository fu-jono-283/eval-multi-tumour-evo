library(lhs)
library(dplyr)

start_time <- Sys.time()
options(warn = 1)

# Defining parameters
params <- data.frame(
  no_of_taxa = c(6, 12),
  individuals_per_taxa = c(4, 8),
  effective_population_size = c(10000, 100000),
  species_birth_rate = c(0.0000001, 0.0001),
  tree_wide_sub_rate = c(0.00000001, 0.00001),
  no_of_sites = c(10000, 100000),
  exponential_growth_rate = c(0.00001, 0.001),
  ado = c(0, 0.2),
  amplification_error = c(0.001, 0.1),
  sequencing_errors = c(0, 0.05)
)

# Mapping the names in R to the parameters for SIMPHY
param_mapping_simphy <- list(
  "no_of_taxa" = "-sl f:%d",
  "individuals_per_taxa" = "-si f:%d",
  "effective_population_size" = "-sp f:%d",
  "species_birth_rate" = "-sb f:%f",
  "tree_wide_sub_rate" = "-su f:%f"
)

# I'm doing the same here for CellCoal
param_mapping_cellcoal <- list(
  "no_of_sites" = "l%d",
  "exponential_growth_rate" = "g%.4e",
  "ado" = "D%.f",
  "amplification_error" = "A%.f",
  "sequencing_errors" = "E%.f",
  "no_of_cells" = "s%d"
)

checkLatinHypercube <- function(X) {
  if (any(apply(X, 2, min) <= 0))
    return(FALSE)
  if (any(apply(X, 2, max) >= 1))
    return(FALSE)
  if (any(is.na(X)))
    return(FALSE)
  # check that the matrix is a Latin hypercube
  g <- function(Y) {
    # check that this column contains all the cells
    breakpoints <- seq(0, 1, length = length(Y) + 1)
    h <- hist(Y, plot = FALSE, breaks = breakpoints)
    all(h$counts == 1)
  }
  # check all the columns
  return(all(apply(X, 2, g)))
}

# Latin Hypercube Sampling design------------------------------------------------------------------------
n <- 2
replicates <- 1

for (i in seq_along(replicates)) {
  initial_lhs <- randomLHS(n, ncol(params))

  total_rows <- nrow(initial_lhs) * replicates
  replicated_lhs <- do.call(rbind, replicate(replicates, as.data.frame(initial_lhs), simplify = FALSE))
    
  # Scaling    
  scaled_lhs <- t(apply(replicated_lhs, 1, function(x) {
    param_range <- unlist(params)
    scaled_vals <- x * (param_range[seq(2, length(param_range), 2)] - param_range[seq(1, length(param_range), 2)]) + param_range[seq(1, length(param_range), 2)]
    scaled_vals[c(1, 2, 3, 6)] <- round(scaled_vals[c(1, 2, 3, 6)])
    scaled_vals
  }))
    
  if (!checkLatinHypercube(initial_lhs)) {
    stop("Generated matrix does not satisfy LHS properties")
  } else {
    cat("SUCCESS: Generated matrix satisfies LHS properties\n")
  }    
    
  colnames(scaled_lhs) <- colnames(initial_lhs)
  param_names <- colnames(params)
  colnames(scaled_lhs) <- param_names
  replicate_labels <- rep(paste0("r", 1:replicates), each = nrow(initial_lhs))

  combined_df <- data.frame(
    scaled_lhs,
    Replicate = as.character(replicate_labels[seq_len(nrow(scaled_lhs))]),
    stringsAsFactors = FALSE
  )

  # Adding score columns 
  scSEQ_score <- rep(-1, nrow(combined_df))
  hapconsensus_score <- rep(-1, nrow(combined_df))
  vcfpseudobulk_score <- rep(-1, nrow(combined_df))
#  allcount_score <- rep(-1, nrow(combined_df))

  combined_df$scSEQ_score <- scSEQ_score
  combined_df$hapconsensus_score <- hapconsensus_score
  combined_df$vcfpseudobulk_score <- vcfpseudobulk_score
#  combined_df$allcount_score <- allcount_score
    
  csv_filename <- sprintf("dataframes/part1_template/combined_data.csv")
  write.csv(combined_df, csv_filename, row.names = FALSE)
}

simphy_folder <- "data/1_simphy"

param_file <- "dataframes/part1_template/combined_data.csv"
param_data <- read.csv(param_file)

# Looping through simphy
for (row in 1:nrow(param_data)) {
  p <- param_data[row, ]
  iteration_identifier <- paste("_pop", p$effective_population_size, "_sites", p$no_of_sites, sep = "")
  no_of_cells <- as.integer(p$individuals_per_taxa) * as.integer(p$no_of_taxa)
  param_mapping_cellcoal[["no_of_cells"]] <- sprintf("s%d", no_of_cells)

  simphy_folder_iteration <- file.path(simphy_folder, iteration_identifier)
  simphy_folder_iteration <- simphy_folder_iteration[!grepl("\\.ipynb_checkpoints", simphy_folder_iteration)]

  config_file <- readLines("scripts/parameter_files/premade_simphy.conf", warn = FALSE)
  config_file <- gsub("-sl f:%d", sprintf("-sl f:%d", p$no_of_taxa), config_file)
  config_file <- gsub("-si f:%d", sprintf("-si f:%d", p$individuals_per_taxa), config_file)
  config_file <- gsub("-sp f:%d", sprintf("-sp f:%d", p$effective_population_size), config_file)
  config_file <- gsub("-sb f:%f", sprintf("-sb f:%.10f", as.numeric(p$species_birth_rate)), config_file)
  config_file <- gsub("-su f:%f", sprintf("-su f:%.10f", as.numeric(p$tree_wide_sub_rate)), config_file)

  config_file <- gsub("-o .*", sprintf("-o %s", file.path(simphy_folder, iteration_identifier)), config_file)
  config_file <- gsub("output_sample", paste0("output_sample_", replicate_labels[row]), config_file)

  updated_config_file <- tempfile()
  writeLines(config_file, updated_config_file)

  cmd_simphy <- paste("simphy -i", updated_config_file)
  system(cmd_simphy)
  file.remove(updated_config_file)

  for (i in 1:8) {
    tree_identifier <- paste("r", i, sep = "")
    source_tree_file <- file.path(simphy_folder_iteration, "1", paste0("g_trees", i, ".trees"))
    destination_tree_folder <- file.path(simphy_folder_iteration, paste0(iteration_identifier, tree_identifier))

    dir.create(destination_tree_folder, showWarnings = FALSE, recursive = TRUE)
    file.copy(source_tree_file, destination_tree_folder)
    file.remove(source_tree_file)

    source_stree_file <- file.path(simphy_folder_iteration, "1", "s_tree.trees")
    destination_stree_file <- file.path(destination_tree_folder, "s_tree.trees")
    file.copy(source_stree_file, destination_stree_file, overwrite = TRUE)

    source_params_file <- file.path(simphy_folder_iteration, paste0(iteration_identifier, ".params"))
    destination_params_file <- file.path(destination_tree_folder, paste0(iteration_identifier, ".params"))
    file.copy(source_params_file, destination_params_file, overwrite = TRUE)

    source_conf_file <- file.path(simphy_folder_iteration, paste0(iteration_identifier, ".conf"))
    destination_conf_file <- file.path(destination_tree_folder, paste0(iteration_identifier, ".conf"))
    file.copy(source_conf_file, destination_conf_file, overwrite = TRUE)

    source_command_file <- file.path(simphy_folder_iteration, paste0(iteration_identifier, ".command"))
    destination_command_file <- file.path(destination_tree_folder, paste0(iteration_identifier, ".command"))
    file.copy(source_command_file, destination_command_file, overwrite = TRUE)
  }

  for (i in 1:8) {
    tree_identifier <- paste("r", i, sep = "")
    destination_tree_folder <- file.path(simphy_folder_iteration, paste0(iteration_identifier, tree_identifier))
    destination_tree_folder_new <- file.path(simphy_folder, paste0(iteration_identifier, tree_identifier))
    file.rename(destination_tree_folder, destination_tree_folder_new)
  }

  if (dir.exists(simphy_folder_iteration)) {
    unlink(simphy_folder_iteration, recursive = TRUE)
  }
}

csv_filenames <- c("dataframes/part1_template/combined_data.csv")

for (csv_filename in csv_filenames) {
  df <- read.csv(csv_filename, stringsAsFactors = FALSE)
  new_df <- data.frame()

  for (replicate_label in unique(df$Replicate)) {
    replicate_df <- df[df$Replicate == replicate_label, ]
    replicated_rows <- lapply(1:8, function(i) {
      replicate_df$Replicate <- paste0("r", i)
      return(replicate_df)
    })

    new_replicate_df <- bind_rows(replicated_rows)
    new_df <- bind_rows(new_df, new_replicate_df)
  }

  write.csv(new_df, csv_filename, row.names = FALSE)
  cat(paste("Updated", csv_filename, "\n"))
}

end_time <- Sys.time()

total_runtime <- end_time - start_time
print(paste("Total runtime:", total_runtime))

# CellCoal------------------------------------------------------------------------
total_jobs <- 16 # adjust based on the amount of required jobs 
max_jobs_per_submission <- 999
num_submissions <- ceiling(total_jobs / max_jobs_per_submission)

template_script_path <- "scripts/2_ccoal_scripts/cellcoal.sh"

for (submission_index in 1:num_submissions) {
  start_index <- (submission_index - 1) * max_jobs_per_submission + 1
  end_index <- min(submission_index * max_jobs_per_submission, total_jobs)

  slurm_script <- readLines(template_script_path)
  slurm_script <- gsub("#SBATCH --array=.*", sprintf("#SBATCH --array=%d-%d", start_index, end_index), slurm_script)
 
  temp_script_file <- tempfile(fileext = ".sh")
  writeLines(slurm_script, con = temp_script_file)

  system(sprintf("sbatch %s", temp_script_file), intern = TRUE, ignore.stderr = TRUE)
  

while (TRUE) {
  current_user <- Sys.getenv("USER")
  squeue_command <- c("squeue", "-u", current_user, "--name=cellcoal_array")
  squeue_output <- system2(squeue_command, stdout = TRUE, stderr = TRUE)
  print(squeue_output)

  if (length(squeue_output) == 0) {
    break
  }

  # Check if there are any jobs in PENDING or RUNNING state
  if (any(grepl("PENDING|RUNNING", squeue_output))) {
    Sys.sleep(60)
  } else {
    break
  }
}
}

# Copying a dataframe to each iteration-------------------------------------------------------------------
source_file <- "dataframes/part1_template/combined_data.csv"

target_directory <- "data/2_ccoal"
subdirectories <- list.dirs(target_directory, full.names = TRUE, recursive = FALSE)

for (subdir in subdirectories) {
  target_file <- file.path(subdir, "combined_data.csv")
  if (!file.exists(target_file)) {
    file.copy(source_file, target_file, overwrite = TRUE)
  
    cat("Copied", source_file, "to", target_file, "\n")
  } else {
    cat("Skipped copying", source_file, "to", target_file, "because it already exists.\n")
  }
}

# Check HCS process ---------------------------------------------------------------------------------
while (sum(sapply(subdirectories, function(subdir) {
  target_path <- file.path(subdir, "full_haplotypes_dir")
  target_file <- file.path(target_path, "full_hap.0001_consensus.fasta")
cat("Contents of the directory:", list.files(target_path), "\n")
  
  file.exists(target_file)
})) != 16) {   #adjust!!!
  counts <- sapply(subdirectories, function(subdir) {
    target_path <- file.path(subdir, "full_haplotypes_dir")
    target_file <- file.path(target_path, "full_hap.0001_consensus.fasta")
  
    file.exists(target_file)
  })
  
  cat("Counts of 'full_hap.0001_consensus.fasta' files in each subdirectory:", counts, "\n")
  
  cat("Waiting for a total of 16 'full_hap.0001_consensus.fasta' files in all subdirectories. Retrying in 1 minute...\n") #adjust!!!
  Sys.sleep(20)
}

# CellPhy Haplotype Consensus------------------------------------------------------------------------
check_condition <- function(bestTree_files) {
  return(length(bestTree_files) == 16) # adjust based on the amount of expected reconstructed trees (same as the number of jobs)
}

script_folder <- "scripts/3_cellphy_scripts/hap_consensus"
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
  Sys.sleep(150)
}

# CellPhy VCF Pseodobulk------------------------------------------------------------------------
check_condition <- function(bestTree_files) {
  return(length(bestTree_files) == 16) # adjust
}

script_folder <- "scripts/3_cellphy_scripts/vcf_consensus"
script_pattern <- "^pre.*\\.R$"
script_files <- list.files(path = script_folder, pattern = script_pattern, full.names = TRUE)

while (TRUE) {
  for (script_file in script_files) {
    cat("Executing script:", script_file, "\n")
    system(paste("Rscript", shQuote(script_file)))
  }

  files_vcf <- list.files("data/2_ccoal", pattern = "^vcf\\.raxml\\.bestTree$", full.names = TRUE, recursive = TRUE)

  cat("List of best tree files:", files_vcf, "\n")

  condition_met <- check_condition(files)
  if (condition_met) {
    cat("All trees have been produced. Exiting the loop.\n")
    break
  }
  Sys.sleep(2000)
}

# CellPhy scSEQ------------------------------------------------------------------------
check_condition <- function(bestTree_files) {
  return(length(bestTree_files) == 16) # adjust 
}

script_folder <- "scripts/3_cellphy_scripts/single_cseq"
script_pattern <- "^pre.*\\.R$"
script_files <- list.files(path = script_folder, pattern = script_pattern, full.names = TRUE)

while (TRUE) {
  for (script_file in script_files) {
    cat("Executing script:", script_file, "\n")
    system(paste("Rscript", shQuote(script_file)))
  }

  files_scseq <- list.files("data/2_ccoal", pattern = "^single\\.raxml\\.bestTree$", full.names = TRUE, recursive = TRUE)

  cat("List of best tree files:", files_scseq, "\n")

  condition_met <- check_condition(files)
  if (condition_met) {
    cat("All trees have been produced. Exiting the loop.\n")
    break
  }
  Sys.sleep(300)
}

# NCOUNT -------------------------------------------------------------------------------------- 



# Compiling the scores -------------------------------------------------------------------------
subdirs <- list.dirs("data/2_ccoal", full.names = TRUE, recursive = FALSE)

template_df <- read.csv("dataframes/part1_template/combined_data.csv")

compiled_hapcon_seq_scores <- rep(-1, nrow(template_df))
compiled_vcfbulk_seq_scores <- rep(-1, nrow(template_df))
compiled_single_seq_scores <- rep(-1, nrow(template_df))
#compiled_allcount_score <- rep(-1, nrow(template_df))

for (subdir in subdirs) {
  csv_file <- file.path(subdir, "combined_data.csv")
  
  if (file.exists(csv_file)) {
    df <- read.csv(csv_file)
    
    updated_hapcon_seq_scores <- df$hapconsensus_score
    updated_vcfbulk_seq_scores <- df$vcfpseudobulk_score
    updated_single_seq_scores <- df$scSEQ_score
      
#   updated_allcount_scores <- df$allcount_score 
    
    updated_hapcon_indices <- which(updated_hapcon_seq_scores >= 0)
    updated_vcfbulk_indices <- which(updated_vcfbulk_seq_scores >= 0)
    updated_single_indices <- which(updated_single_seq_scores >= 0)
      
#    updated_allcount_indices <- which(updated_allcount_scores >= 0)

    compiled_hapcon_seq_scores[updated_hapcon_indices] <- updated_hapcon_seq_scores[updated_hapcon_indices]
    compiled_vcfbulk_seq_scores[updated_vcfbulk_indices] <- updated_vcfbulk_seq_scores[updated_vcfbulk_indices]
    compiled_single_seq_scores[updated_single_indices] <- updated_single_seq_scores[updated_single_indices]
      
#    compiled_allcount_scores[updated_allcount_indices] <- updated_allcount_scores[updated_allcount_indices]  
  }
}

template_df$vcfpseudobulk_score <- compiled_vcfbulk_seq_scores
template_df$scSEQ_score <- compiled_single_seq_scores
template_df$hapconsensus_score <- compiled_hapcon_seq_scores

#template_df$allcount_score <- compiled_allcount_scores

write.csv(template_df, "dataframes/testing_df/part_one.csv", row.names = FALSE)
