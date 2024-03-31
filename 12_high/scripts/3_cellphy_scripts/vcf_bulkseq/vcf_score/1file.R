initial_path <- "data/2ccoal"
iteration_identifier <- "_pop99731_sites74578r8"

vcf_dir_path <- file.path(initial_path, iteration_identifier, "vcf_dir")

if (any(grepl("vcf.raxml.bestTree", list.files(file.path(initial_path, iteration_identifier))))) {
  cat("Skipping iteration:", iteration_identifier, " - folder contains vcf.raxml files\n")
} else if (!file.exists(file.path(initial_path, iteration_identifier, "vcf_dir", "vcf_no_outgcell.recode.vcf_consensus.vcf"))) {
  cat("Skipping iteration:", iteration_identifier, " - folder does not contain the required file\n")
} else {
  num_sites <- as.numeric(gsub(".*_sites(\\d+).*", "\\1", iteration_identifier))

  if (num_sites > 50000) {
    cpus_per_task <- 1
    mem <- "8G"  
  } else {
    cpus_per_task <- 1
    mem <- "8G"  
  }

  cellphy_vcf_sl <- readLines("scripts/3_cellphy_scripts/vcf_bulkseq/vcf_score/vcf_score.sh")

  modified_vcf_script <- gsub("--cpus-per-task=\\d+", paste0("--cpus-per-task=", cpus_per_task), cellphy_vcf_sl)
  modified_vcf_script <- gsub("--mem=\\d+G", paste0("--mem=", mem), modified_vcf_script)  
  modified_vcf_script <- gsub("\\$cpus_per_task", cpus_per_task, modified_vcf_script)
  modified_vcf_script <- gsub("\\$iteration_identifier", iteration_identifier, modified_vcf_script)
    
  writeLines(modified_vcf_script, "scripts/3_cellphy_scripts/vcf_bulkseq/vcf_score/mod_vcf_score.sh")

  cmd_cellphy_vcf <- paste("sbatch scripts/3_cellphy_scripts/vcf_bulkseq/vcf_score/mod_vcf_score.sh", sep="")
      
  system(cmd_cellphy_vcf)
}
