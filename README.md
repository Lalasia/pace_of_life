# The-pace-of-life-for-forest-trees
Code for the paper titled: The pace of life for forest trees

Bialic-Murphy, et al 2024

--------------------------------------------------------------
Contents of "script" folder:
--------------------------------------------------------------
1_Fit_Vital_R  --Load data and fit vital rate functions.
2_IPM_4_HPC.R  --Load vital rate parameter estimates, set up IPMs, calculate LHTs 
2_IPM_4_HPC.sh  --Sets up 2_IPM_4_HPC.R to run on an HPC.
3_IPM_postProc.R  --Compiles all runs generated from 2_IPM_4_HPC.R. Creates data files for subsequent scripts
4_Clusters_analysis.R  --Runs cluster analysis and generates figures
5_Bayesian_model.R  --Runs analysis of environmental effects on LHTs
6_Hull_volume.R --Calculates the convex hull of LHTs and analyses patterns
Error_Checks --Contains modified versions of 2_IPM_4_HPC.sh and 2_IPM_4_HPC.R (See Note 2 below).

--Run these scripts to recreate our analyses and results.
--The input data required for these scripts and the intermediate output data that are both produced and later required by these scripts are stored on Zanodo: doi.org/10.5281/zenodo.11615767

--NOTE 1 | The raw tree-by-tree obs. data from networks with sensitive species are not included here but can be accessed by submitting a data request to the associated networks. See our DMA statement for more details. Here, we provide the raw tree-by-tree data from the fully open-access networks. This affects the first R script only. The grid-level output data, generated in script 1 and shared in the data folder, can be used to fully replicate the results of this study.
	
--NOTE 2 | The bash script [2_IPM_4_HPC.sh] sets up an array to run the IPM script [2_IPM_4_HPC.R] in batches of 100 species on a high-performance computing cluster.
Note that two species [index=30 & 962] had an immortal survival-growth matrix. To find this error, we duplicated the IPM R and bash scripts and modified the indexing. This can be replicated using the IPM_4_Euler R scripts in the "Error Check" folder.
	
3_IPM_postProc.R is then used to recompile the IPM outputs.

--------------------------------------------------------------
Contents of "data" folder:
--------------------------------------------------------------
Initial input data are stored on Zanodo: doi.org/10.5281/zenodo.11615767

TreeMort_TreeData.Rdata --Cannot share original data due to third-party data restrictions. The included subset are the data from networks that are fully open access. The other datasets can be accessed by submitting a data request to the networks. See DMA statement for more details
hex_250k.csv --Contains climate data at grid cell level. Data coordinates are not provided due to sensitive species.

Outputs from Script 1 and grid-level input data needed to replicate models and analysis, i.e., Script 2 - 6
TreeMort_Growth_1yr_woGrid.Rdata --Growth function parameters
TreeMort_sigmag_1yr_woGrid.Rdata --Variance estimates for growth function
TreeMort_Survival_1yr_woGrid.Rdata --Survival function parameters
Min_Max_Size.Rdata --Species size characteristics
