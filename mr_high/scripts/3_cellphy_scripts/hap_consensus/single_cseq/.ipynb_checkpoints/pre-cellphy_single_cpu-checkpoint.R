# The script for single_cseq is slightly different in that it requires a term called result_variable to be considered instead of just the number of sites. The result variable is the product of the number of sites multiplied by the number of individuals... Therefore the total number of sites in the full haplotype file is considered when adjusting for the amount of cpu to be used for the Cellphy analysis. 

iteration_identifiers <- system(paste("ls -1", shQuote("data/2_ccoal"), "| sed -n", Sys.getenv("SLURM_ARRAY_TASK_ID"), "p"), intern = TRUE)

for (iteration_identifier in iteration_identifiers) {
  if (any(grepl("^single.raxml.bestTree", list.files(file.path("data/2_ccoal", iteration_identifier))))) {
    cat("Skipping iteration:", iteration_identifier, " - folder contains single.raxml files\n")
    next
  }

  full_hap_file <- file.path("data/2_ccoal", iteration_identifier, "full_haplotypes_dir", "full_hap.0001")

  if (file.exists(full_hap_file)) {
    full_hap_data <- readLines(full_hap_file, n = 1)
    num_from_full_hap <- as.numeric(sub("^[^0-9]*([0-9]+).*", "\\1", full_hap_data))
    cat("First number from full_hap.0001:", num_from_full_hap, "\n")
  } else {
    cat("full_hap.0001 not found in iteration:", iteration_identifier, "\n")
  }

  num_sites <- as.numeric(gsub(".*_sites(\\d+).*", "\\1", iteration_identifier))

  result_variable <- num_sites * num_from_full_hap

  if (result_variable > 30000000) {
    cat("Skipping iteration:", iteration_identifier, " - result_variable is greater than set limit\n")
    next
  }

  if (result_variable > 20000000) {
    cpus_per_task <- 4
    mem <- "20G"  
  } else if (result_variable > 10000000) {
    cpus_per_task <- 4
    mem <- "20G"    
   } else if (result_variable > 7500000) {
    cpus_per_task <- 4
    mem <- "10G"     
  } else if (result_variable > 5000000) {
    cpus_per_task <- 4
    mem <- "10G"  
  } else if (result_variable > 2500000) {
    cpus_per_task <- 4
    mem <- "10G"
  } else if (result_variable > 1000000) {
    cpus_per_task <- 4
    mem <- "10G" 
  } else if (result_variable > 500000) {
    cpus_per_task <- 4
    mem <- "4G" 
  } else {
    cpus_per_task <- 4
    mem <- "4G"  
  }

  cellphy_single_sh <- readLines("scripts/3_cellphy_scripts/single_cseq/cellphy_single.sh")

  modified_single_script <- gsub("--cpus-per-task=\\d+", paste0("--cpus-per-task=", cpus_per_task), cellphy_single_sh)
  modified_single_script <- gsub("--mem=\\d+G", paste0("--mem=", mem), modified_single_script) 
  modified_single_script <- gsub("\\$cpus_per_task", cpus_per_task, modified_single_script)
  modified_single_script <- gsub("\\$iteration_identifier", iteration_identifier, modified_single_script)

  writeLines(modified_single_script, "scripts/3_cellphy_scripts/single_cseq/mod_cellphy_single.sh")

  cmd_cellphy_single <- paste("sbatch scripts/3_cellphy_scripts/single_cseq/mod_cellphy_single.sh", sep="")

  system(cmd_cellphy_single)
}
