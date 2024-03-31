iteration_identifiers <- system(paste("ls -1", shQuote("data/2_ccoal"), "| sed -n", Sys.getenv("SLURM_ARRAY_TASK_ID"), "p"), intern = TRUE)

for (iteration_identifier in iteration_identifiers) {
  if (any(grepl("^full_hap_consensus_1.00.fasta", list.files(file.path("data/2_ccoal", iteration_identifier, "full_haplotypes_dir"))))) {
    cat("Skipping iteration:", iteration_identifier, " - folder contains the last cutoff consensus_seq file\n")
    next
  }
    
  num_sites <- as.numeric(gsub(".*_sites(\\d+).*", "\\1", iteration_identifier))

  if (num_sites > 300000) {
    cpus_per_task <- 2
    mem <- "4G"  
  } else if (num_sites > 250000) {
    cpus_per_task <- 2
    mem <- "4G"  
  } else if (num_sites > 100000) {
    cpus_per_task <- 2
    mem <- "4G" 
  } else if (num_sites > 50000) {
    cpus_per_task <- 2
    mem <- "4G" 
  } else {
    cpus_per_task <- 2
    mem <- "4G"  
  }
    
cutoff_sh <- readLines("scripts/count/cutoff/cutoff.sh")

modified_cutoff_script <- gsub("--cpus-per-task=\\d+", paste0("--cpus-per-task=", cpus_per_task), cutoff_sh)
modified_cutoff_script <- gsub("--mem=\\d+G", paste0("--mem=", mem), modified_cutoff_script)  
modified_cutoff_script <- gsub("\\$cpus_per_task", cpus_per_task, modified_cutoff_script)

writeLines(modified_cutoff_script, "scripts/count/cutoff/mod_cutoff.sh")

cmd_cutoff <- paste("sbatch --export=iteration_identifier=\"", iteration_identifier, "\" scripts/count/cutoff/mod_cutoff.sh", sep="")

system(cmd_cutoff)
}