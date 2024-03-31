# For re-processing the consensus.fasta files after adjusting for optimal cutoff value. 

iteration_identifiers <- system(paste("ls -1", shQuote("data/2_ccoal"), "| sed -n", Sys.getenv("SLURM_ARRAY_TASK_ID"), "p"), intern = TRUE)

for (iteration_identifier in iteration_identifiers) {
  if (any(grepl("^full_hap.0001_consensus.fasta", list.files(file.path("data/2_ccoal", iteration_identifier, "full_haplotypes_dir"))))) {
    cat("Skipping iteration:", iteration_identifier, " - folder contains the consensus.fasta file already\n")
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
    
consensus_redo <- readLines("scripts/2_ccoal_scripts/consensus_seq_redo/consensus_redo.sh")

modified_con_script <- gsub("--cpus-per-task=\\d+", paste0("--cpus-per-task=", cpus_per_task), consensus_redo)
modified_con_script <- gsub("--mem=\\d+G", paste0("--mem=", mem), modified_con_script)  
modified_con_script <- gsub("\\$cpus_per_task", cpus_per_task, modified_con_script)

writeLines(modified_con_script, "scripts/2_ccoal_scripts/consensus_seq_redo/mod_consensus_redo.sh")

cmd_con <- paste("sbatch --export=iteration_identifier=\"", iteration_identifier, "\" scripts/2_ccoal_scripts/consensus_seq_redo/mod_consensus_redo.sh", sep="")

system(cmd_con)
}
