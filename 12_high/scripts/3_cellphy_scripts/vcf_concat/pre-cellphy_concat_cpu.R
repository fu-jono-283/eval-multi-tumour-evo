iteration_identifiers <- system(paste("ls -1", shQuote("data/7_vcf_concat"), "| sed -n", Sys.getenv("SLURM_ARRAY_TASK_ID"), "p"), intern = TRUE)

for (iteration_identifier in iteration_identifiers) {
  if (any(grepl("^vcf_concat.raxml.bestTree", list.files(file.path("data/7_vcf_concat", iteration_identifier))))) {
    cat("Skipping iteration:", iteration_identifier, " - folder contains vcf_concat.raxml files\n")
    next
  }

  num_sites <- as.numeric(gsub(".*_sites(\\d+).*", "\\1", iteration_identifier))

  if (num_sites > 300000) {
    cpus_per_task <- 10
    mem <- "20G"  
  } else if (num_sites > 250000) {
    cpus_per_task <- 10
    mem <- "20G"  
  } else if (num_sites > 100000) {
    cpus_per_task <- 10
    mem <- "20G" 
  } else if (num_sites > 50000) {
    cpus_per_task <- 10
    mem <- "10G" 
  } else {
    cpus_per_task <- 10
    mem <- "10G"  
  }

cellphy_concat_sh <- readLines("scripts/3_cellphy_scripts/vcf_concat/cellphy_concat.sh")

modified_concat_script <- gsub("--cpus-per-task=\\d+", paste0("--cpus-per-task=", cpus_per_task), cellphy_concat_sh)
modified_concat_script <- gsub("--mem=\\d+G", paste0("--mem=", mem), modified_concat_script)     
modified_concat_script <- gsub("\\$cpus_per_task", cpus_per_task, modified_concat_script)
modified_concat_script <- gsub("\\$iteration_identifier", iteration_identifier, modified_concat_script)

writeLines(modified_concat_script, "scripts/3_cellphy_scripts/vcf_concat/mod_cellphy_concat.sh")

cmd_cellphy_concat <- paste("sbatch scripts/3_cellphy_scripts/vcf_concat/mod_cellphy_concat.sh", sep="")

system(cmd_cellphy_concat)
}
