# R packages
library(lhs)
library(dplyr)


# Recording duration of simulation 
start_time <- Sys.time()
options(warn = 1)

# Defining parameters of LHS dataframe
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

# Mapping the names in R to the parameters for SIMPHY software
param_mapping_simphy <- list(
  "no_of_taxa" = "-sl f:%d",
  "individuals_per_taxa" = "-si f:%d",
  "effective_population_size" = "-sp f:%d",
  "species_birth_rate" = "-sb f:%f",
  "tree_wide_sub_rate" = "-su f:%f"
)

# Doing the same here for CellCoal software
param_mapping_cellcoal <- list(
  "no_of_sites" = "l%d",
  "exponential_growth_rate" = "g%.4e",
  "ado" = "D%.f",
  "amplification_error" = "A%.f",
  "sequencing_errors" = "E%.f",
  "no_of_cells" = "s%d"
)

# Function that ensures that LHS is being generated correctly
checkLatinHypercube <- function(X) {
  if (any(apply(X, 2, min) <= 0))
    return(FALSE)
  if (any(apply(X, 2, max) >= 1))
    return(FALSE)
  if (any(is.na(X)))
    return(FALSE)
  g <- function(Y) {
    breakpoints <- seq(0, 1, length = length(Y) + 1)
    h <- hist(Y, plot = FALSE, breaks = breakpoints)
    all(h$counts == 1)
  }
  return(all(apply(X, 2, g)))
}

# Latin Hypercube Sampling design---------------------------------------------------------

n <- 200 # generating 200 samples  
replicates <- 1 # not performing augment LHS - just one full simulation of dataframe here

# Forming dataframe
for (i in seq_along(replicates)) {
  initial_lhs <- randomLHS(n, ncol(params))

  total_rows <- nrow(initial_lhs) * replicates
  replicated_lhs <- do.call(rbind, replicate(replicates, as.data.frame(initial_lhs), simplify = FALSE))
    
# Scaling of LHS values from 0-1 to those of actual parameters   
  scaled_lhs <- t(apply(replicated_lhs, 1, function(x) {
    param_range <- unlist(params)
    scaled_vals <- x * (param_range[seq(2, length(param_range), 2)] - param_range[seq(1, length(param_range), 2)]) + param_range[seq(1, length(param_range), 2)]
    scaled_vals[c(1, 2, 3, 6)] <- round(scaled_vals[c(1, 2, 3, 6)])
    scaled_vals
  }))

# Implementing the LHS check function    
  if (!checkLatinHypercube(initial_lhs)) {
    stop("Generated matrix does not satisfy LHS properties")
  } else {
    cat("SUCCESS: Generated matrix satisfies LHS properties\n")
  }    

# Naming columns of dataframe, adding in replicate labels    
  colnames(scaled_lhs) <- colnames(initial_lhs)
  param_names <- colnames(params)
  colnames(scaled_lhs) <- param_names
  replicate_labels <- rep(paste0("r", 1:replicates), each = nrow(initial_lhs))

# Finalizing as combined_df, which includes the concept of replicates within the initial df. 
  combined_df <- data.frame(
    scaled_lhs,
    Replicate = as.character(replicate_labels[seq_len(nrow(scaled_lhs))]),
    stringsAsFactors = FALSE
  )

# Adding score columns 
  scSEQ_score <- rep(-1, nrow(combined_df))
  hapcon_seq_score <- rep(-1, nrow(combined_df))
  vcfbulk_seq_score <- rep(-1, nrow(combined_df))

  combined_df$scSEQ_score <- scSEQ_score
  combined_df$hapcon_seq_score <- hapcon_seq_score
  combined_df$vcfbulk_seq_score <- vcfbulk_seq_score
  
# Writing dataframe into a CSV file - based on number of samples
  csv_filename <- sprintf("dataframes/part1_template/combined_data_%02d.csv", n)
  write.csv(combined_df, csv_filename, row.names = FALSE)
}

# Specifying where input/output of SimPhy data will go
simphy_folder <- "data/1_simphy"

param_file <- "dataframes/part1_template/combined_data_200.csv" # Needs adjustment based on number of samples simulated
param_data <- read.csv(param_file)

# Looping through simphy
for (row in 1:nrow(param_data)) {
  p <- param_data[row, ]
# The iteration identifier - a concept used throughout simulations to identify each iteration - based initially
# on the EPS, NOS and later based also on the replicate number (r1-r8)
  iteration_identifier <- paste("_pop", p$effective_population_size, "_sites", p$no_of_sites, sep = "")
# The number of cells is a parameter argument present in CellCoal - which must be calculated from the SimPhy values
# simulated (individuals per taxa, and number of taxa).
  no_of_cells <- as.integer(p$individuals_per_taxa) * as.integer(p$no_of_taxa)
  param_mapping_cellcoal[["no_of_cells"]] <- sprintf("s%d", no_of_cells)

  simphy_folder_iteration <- file.path(simphy_folder, iteration_identifier)
  simphy_folder_iteration <- simphy_folder_iteration[!grepl("\\.ipynb_checkpoints", simphy_folder_iteration)]
# Reading the OG configuration file of SimPhy, contains the execution script made by SimPhy creators, in which we
# susbtitute our LHS values
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

# Note that in SimPhy we generate 8 gene trees per sample... Now we specify the iteration identifier to consider each replicate 
# as an individual iteration. In these loops, we rearrange our folders/files so that each replicate has its own output folder, 
# and each output folder contains the gene tree, species tree etc. 
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

# We update the dataframe so that it considers that there are now 8 replicates - so the number of rows equates to 
# no. of samples * 8. 
csv_filenames <- c("dataframes/part1_template/combined_data_200.csv")

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

# Recording simulation time
end_time <- Sys.time()

total_runtime <- end_time - start_time
print(paste("Total runtime:", total_runtime))
