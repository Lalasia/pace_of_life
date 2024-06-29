#!/bin/bash
# Job script to run an R script on a single node
#
#SBATCH --job-name=IPM_trial_1         # Job Name
#SBATCH --mem-per-cpu=24G              #request 12GB of RAM
#SBATCH --qos normal                   # normal priority level
#SBATCH --ntasks=1                     # Run on a single core
#SBATCH --time=2-05:00                 # Time limit days-hrs:min
#SBATCH -o ./Output/IPM_%A_%a.out # Standard output - write the %A= Job ID and %a = task or Step ID
#SBATCH -e ./Output/IPM_%A_%a.err # Standard error 
#SBATCH --array=1-12%12 # create a array for each batch to run SLURM_ARRAY_TASK_ID
#
#load all necessary modules:
#module purge   # clear any inherited modules
module load  r/4.1.3   # load some modules e.g  a specific R version.

#execute your script:
Rscript --vanilla ./Rscript/2_IPM_4_HPC.R ${SLURM_ARRAY_TASK_ID} #Run my R script for given year

#Script requires the following files loaded to HPC account:
#./Rdata/TreeMort_Survival_1yr_woGrid.Rdata
#./Rdata/TreeMort_Growth_1yr_woGrid.Rdata
#./Rdata/TreeMort_sigmag_1yr_woGrid.Rdata
#./Rdata/Min_Max_Size.Rdata