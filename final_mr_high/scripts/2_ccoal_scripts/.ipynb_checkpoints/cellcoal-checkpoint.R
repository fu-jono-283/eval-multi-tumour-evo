# Executed by cellcoal.sh through array jobs, this script processes the CellCoal samples -- removes the outgroup, and also creates a consensus sequence fasta file using a specified cutoff value. The number of samples processed can be adjusted by changing the dataframe that the script initially reads.

start_time <- Sys.time()
options(warn=1)

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

# CELLCOAL---------------------------------------------------------
param_file_cellcoal <- readLines("scripts/configuration_folder/premade_ccoal.param", warn = FALSE)

# Change the combined_data_[ ].csv---------------------------------
combined_df <- read.csv("dataframes/part1_template/combined_data.csv", header = TRUE)
array_job_index <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
p <- combined_df[array_job_index, ]

  effective_population_size <- p$effective_population_size 
         tree_wide_sub_rate <- p$tree_wide_sub_rate
                no_of_sites <- p$no_of_sites
    exponential_growth_rate <- p$exponential_growth_rate
                        ado <- p$ado 
        amplification_error <- p$amplification_error
          sequencing_errors <- p$sequencing_errors
                no_of_cells <- p$no_of_cell
                  Replicate <- p$Replicate

# Where processed iteration folders of CellCoal will go 
single_ccoal_output_folder <- "data/2_ccoal"
iteration_identifier <- paste("_pop", effective_population_size, "_sites", no_of_sites, Replicate, sep = "")

# It will skip the iterations that already have a matching output folder in 2_ccoal_folder
output_folder <- file.path(single_ccoal_output_folder, iteration_identifier)
if (file.exists(output_folder)) {
cat("Skipping iteration:", iteration_identifier, " - Output folder already exists\n")
} else {
    
  no_of_cells <- as.integer(p[2]) * as.integer(p[1])
  param_mapping_cellcoal[["no_of_cells"]] <- sprintf("s%d", no_of_cells)

  param_file_cellcoal <- gsub("e%d", sprintf("e%d", effective_population_size), param_file_cellcoal)
  param_file_cellcoal <- gsub("u%f", sprintf("u%.10f", tree_wide_sub_rate), param_file_cellcoal)
  param_file_cellcoal <- gsub("l%d", sprintf("l%d", no_of_sites), param_file_cellcoal)
  param_file_cellcoal <- gsub("g%.10e", sprintf("g%.10e", exponential_growth_rate), param_file_cellcoal)
  param_file_cellcoal <- gsub("D%f", sprintf("D%.10f", ado), param_file_cellcoal)
  param_file_cellcoal <- gsub("A%f", sprintf("A%.10f", amplification_error), param_file_cellcoal)
  param_file_cellcoal <- gsub("E%f", sprintf("E%.10f", sequencing_errors), param_file_cellcoal)
  param_file_cellcoal <- gsub("s%d", sprintf("s%d", no_of_cells), param_file_cellcoal)

  param_file_cellcoal <- gsub("o%", sprintf("o%s", file.path(single_ccoal_output_folder, iteration_identifier)), param_file_cellcoal)
  param_file_cellcoal <- gsub("T%", sprintf("T%s", file.path(grep("g_trees[0-9]+\\.trees", list.files(file.path("data/1_simphy", iteration_identifier), full.names = TRUE), value = TRUE)[1])), param_file_cellcoal)

temp_param_file_cellcoal <- tempfile(fileext = ".txt")
  writeLines(param_file_cellcoal, temp_param_file_cellcoal)
  cat("Temporary parameter file saved:", temp_param_file_cellcoal, "\n")
  
  cmd_cellcoal <- paste("cellcoal -F", temp_param_file_cellcoal, sep = "")
  system(cmd_cellcoal)
    
  cat("Cellcoal iteration folder processing finished.\n")   
    
# REMOVE THE OUTGCELL---------------------------------------------------------
update_phylip_file <- function(file_path) {
    content <- readLines(file_path)
    outgcell_line <- grep("outgcell", content)
    if (length(outgcell_line) > 0) {
        content <- content[1:(outgcell_line - 1)]
}
    num_taxa <- length(content) - 1
    first_line <- strsplit(content[1], "\\s+")
    first_line[[1]][1] <- as.character(num_taxa)
    content[1] <- paste(first_line[[1]], collapse = " ")
    writeLines(content, file_path)
}

full_hap_file <- file.path(single_ccoal_output_folder, iteration_identifier, "full_haplotypes_dir", "full_hap.0001") 

    content <- readLines(full_hap_file)
    update_phylip_file(full_hap_file)
    print(paste("Outgcell removed from:", full_hap_file))

cat("Outgcell removed.\n")
    
# BUILDING THE CONSENSUS SEQUENCE----------------------------------------------------
cmd_python <- sprintf("python %s %s", "scripts/2_ccoal_scripts/consensus_cutoff.py", full_hap_file)
system(cmd_python)
cat("Consensus sequence file finished.\n")    
}

end_time <- Sys.time()
total_runtime <- end_time - start_time
print(paste("Total runtime:", total_runtime))