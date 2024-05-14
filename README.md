# Evaluating the Use of Bulk and Single-cell Sequences in Phylogenetic Analyses of Multi-tumour Evolution

## Executing simulation

To re-run simulation, please utilise the scripts present in the final_mr_high folder. The two pipelines featured in this folder are used to perform phylogenetic reconstruction using several methods. 

To run part1 requires the installation of SimPhy (Mallo, 2016), CellCoal (Posada, 2020), CellPhy (Kozlov, 2022). 

To run part2 requires the installation of ASTRAL-III (Zhang, 2018), FastRFS (Vachaspati, 2017), CellPhy (Kozlov, 2022).

Adjustments of parameter values can be made in scripts/configuration_folder/ - 1 file for CellCoal, 1 for SimPhy. 

Adjustments of number of samples simulated (number of iteratations and trees must also be changed accordingly) can be made in 1_all_script/1_all.R

To execute the simulation for part 1, load R and enter: Rscript scripts/1_all_script/1_all.R. 

To execute the simulation for part 2... 

Please refer to the [Wiki](https://github.com/fu-jono-283/eval-multi-tumour-evo/wiki/Overview) for an overview of the project. 
