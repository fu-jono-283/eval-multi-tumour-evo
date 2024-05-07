iteration_identifiers <- system(paste("ls -1", shQuote("data/2_ccoal"), "| sed -n", Sys.getenv("SLURM_ARRAY_TASK_ID"), "p"), intern = TRUE)

for (iteration_identifier in iteration_identifiers) {
  if (any(grepl("^hapcon.raxml.bestTree", list.files(file.path("data/2_ccoal", iteration_identifier))))) {
    cat("Skipping iteration:", iteration_identifier, " - folder contains hapcon.raxml files\n")
    next
 }

  num_sites <- as.numeric(gsub(".*_sites(\\d+).*", "\\1", iteration_identifier))

  if (num_sites > 100000) {
    cpus_per_task <- 25
    mem <- "4G"  
  } else if (num_sites > 75000) {
    cpus_per_task <- 15
    mem <- "4G"  
  } else if (num_sites > 50000) {
    cpus_per_task <- 10
    mem <- "4G" 
  } else {
    cpus_per_task <- 6
    mem <- "4G"  
  }

cellphy_hapcon_sh <- readLines("scripts/3_cellphy_scripts/hap_consensus/cellphy_hapcon.sh")

modified_hapcon_script <- gsub("--cpus-per-task=\\d+", paste0("--cpus-per-task=", cpus_per_task), cellphy_hapcon_sh)
modified_hapcon_script <- gsub("--mem=\\d+G", paste0("--mem=", mem), modified_hapcon_script)     
modified_hapcon_script <- gsub("\\$cpus_per_task", cpus_per_task, modified_hapcon_script)
modified_hapcon_script <- gsub("\\$iteration_identifier", iteration_identifier, modified_hapcon_script)

writeLines(modified_hapcon_script, "scripts/3_cellphy_scripts/hap_consensus/mod_cellphy_hapcon.sh")

cmd_cellphy_hapcon <- paste("sbatch scripts/3_cellphy_scripts/hap_consensus/mod_cellphy_hapcon.sh", sep="")

system(cmd_cellphy_hapcon)
}
