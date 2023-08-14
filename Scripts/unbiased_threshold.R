# unbiased_threshold.R: Finds unbiased thresholds for optimal case control distinction

'Usage: 
  unbiased_threshold.R --model=<model> --min_sample=<min_sample> --resolution_var=<resolution_var> --nperm=<nperm> --threshold_var=<threshold_var> [--case_group=<case_group>] [--cores=<cores>]
  
  Options:
  -h --help
  --model=<model> model to use for analysis
  --min_sample=<min_sample> minimum size of case/control 
  --resolution_var=<resolution_var> resolution for clustering.
  --nperm=<nperm> number of permutations
  --threshold_var=<threshold_var> thresholding variable
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
  --cores=<cores> cores to use for permutations. Do not use more than 20 unless other members of group aware. [default: 5]
  
' -> doc
# library(tidyverse)
# library(data.table)
library(here)
library(docopt)
source(here("Scripts/fp_forest_plot_functions.R"))
source(here("Scripts/unbiased_threshold_functions_repo.R"))
arguments <- docopt(doc, version = 'unbiased_threshold.R 2.2')

# debugging
# arguments <- list()
# arguments$model <- "URMIS_min80"
# arguments$min_sample <- "20"
# arguments$resolution_var <- "0_2"
# arguments$nperm <- "100"
# arguments$case_group <- "total_cohort"
# arguments$threshold_var <- "mis_z"
# arguments$cores <- "1"




# Start Log
dir.create(log_folder <- here("Data/unbiased_thresholdLog"))
time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_",arguments$case_group, "_")
logr::log_open(paste0(log_folder,"/", time_case_prefix,"unbiased_threshold.log"))
logr::log_print(arguments)
# Load in genotype file
geno_ped_df <- cluster_sample_and_genotype_list_fp(paste0(here(),"/"),arguments$resolution_var, as.numeric(arguments$min_sample), arguments$model, nclust_overide=NULL,matrix_geno = "geno") #as.numeric(arguments$cores),

if(arguments$threshold_var == "LOEUF"){
  
  dir.create(output_str <- here("Results/loeuf_perms_output/"))
  dir.create(permutation_directory <- here(paste0("Results/loeuf_perms_res_", arguments$resolution, "_min_sample_", arguments$min_sample, "/")))
} else if(arguments$threshold_var == "mis_z"){
  dir.create(output_str <- here("Results/mis_z_perms_output/"))
  dir.create(permutation_directory <- here(paste0("Results/mis_z_perms_res_", arguments$resolution, "_min_sample_", arguments$min_sample, "/")))
}
figure_gene_loeuf_intolerance_threshold_repo(geno_ped_df$genotype_df,geno_ped_df$sample_df,
                                             paste0(arguments$case_group,"_",arguments$model, "_min_case_",arguments$min_sample,"_resolution_",
                                                    arguments$resolution_var),permutation_directory,
                                             output_str,as.numeric(arguments$nperm),as.numeric(arguments$cores),arguments$threshold_var,save_stat = "stat_p", test_flag = "cmh")

# Close log
logr::log_print("Finished")
logr::log_close()


# 




# figure_gene_loeuf_intolerance_threshold(df_list_picu_all_ptv_flex_p9$genotype_df,df_list_picu_all_ptv_flex_p9$sample_df,
#                                         paste0("PICU_all_loeuf_FlexPTV_min_case_pext9_",min_case_picu_all,"_resolution_",resolution_var_picu_all),permutation_directory,
#                                         output_str=output_str,nperm=100000,save_stat = "stat_p", test_flag = "cmh")



