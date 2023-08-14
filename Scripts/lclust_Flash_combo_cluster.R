# lclust_Flash_combo_cluster.R: creates a large combined cluster with multiple ancestries but keeps ratios intact

'Usage: 
  lclust_Flash_combo_cluster.R --resolution_var=<resolution_var> [--min_sample=<min_sample>]  [--case_group=<case_group>] [--max_ratio=<max_ratio>] [--debug=<debug>]
  
  Options:
  -h --help
  --resolution_var=<resolution_var> resolution for clustering. Typically 0.1-0.4
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
  --min_sample=<min_sample> minimum size of case/control allowed to include cluster into combined cluster [default: 0]
  --max_ratio=<max_ratio> desired case/control ratio. Will downsample controls to achieve ratio. If blank, program will find the max value in among clusters and then choose that value. If users sets this, will ignore clusters that are greater than user value and downsample controsl for other clusters [default: FALSE]
  --debug=<debug> [default: FALSE]
  
' -> doc
library(tidyverse)
library(data.table)
# library(stringr)
library(here)
library(docopt)
arguments <- docopt(doc, version = 'lclust_Flash_combo_cluster.R 1.1')

# debugging
if(arguments$debug == "TRUE"){
  arguments <- list()
  arguments$resolution_var <- "0_4"
  arguments$case_group <- "ic_20230523"
  arguments$min_sample <- "100"
  arguments$max_ratio <- "TRUE"
  # arguments$mixed_ancestry_cluster <- "9"
}


# Create log files
tryCatch({
  dir.create(log_folder <- here("Data","lclust_Flash_combo_cluster_log"))
  time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_")
  logr::log_open(here(log_folder,paste0(time_case_prefix,"lclust_Flash_combo_cluster_logfile.log")))
  logr::log_print(arguments)
}, error = function(e) {
  # Handle the error here
  print(paste("An error occurred while creating the logfile:", e$message))
})

# cast variables
resolution <- arguments$resolution_var
min_sample <- as.numeric(arguments$min_sample)

# establish paths
dir <- ""

if (dir == "") {
  dir = ""
} else {
  dir <- paste0(dir, "_")
}

plot_path <- here(paste0("Results/Plots_", gsub("_", "",  dir), "/"))
plot_path <- gsub("_/", "/", plot_path)
results_path <- here(paste0("Results/CMH_", gsub("_", "",  dir), "/"))
results_path <- gsub("_/", "/", results_path)
KinshipFlashPCA_path <- here(paste0("Results/KinshipFlashPCA_", gsub("_", "",  dir), "/"))
KinshipFlashPCA_path <- gsub("_/", "/", KinshipFlashPCA_path)


# load in cluster sizes 2022_02_22_10_41_54_temp_case_lclustering_res_0_2_cluster_sizes_ratios.txt
logr::log_print(cl_size_path <- list.files(plot_path, pattern = paste0("*lclustering_res_", resolution, "_cluster_sizes_ratios.txt$"), full.names = TRUE))
cl_sizes <- fread(cl_size_path)
logr::log_print(cl_sizes)
logr::log_print("Limiting clusters based on acceptable case and control sample sizes")
logr::log_print(nclust_size <- (cl_sizes %>% filter(case >= min_sample & control >= min_sample))$cluster %>% print())
cl_sizes <- cl_sizes %>% mutate(ratio_josh = case/control)

# find max case to control ratio
max_ratio_possible <- max((cl_sizes %>% filter(cluster %in%  nclust_size))$ratio_josh) %>% print()
logr::log_print(paste0("Largest case/control ratio possible: ", max_ratio_possible))
if(arguments$max_ratio == "FALSE"){
  max_ratio <- max_ratio_possible
  logr::log_print(paste0("User did not define max ratio so selecting from clusters with both case"))
} else {
  logr::log_print(paste0("User selected case/control ratio: ", max_ratio <- max_ratio_possible))
  # logr::log_print(paste0("Case/control ratio used for this max cluster: ", max_ratio))
}

logr::log_print("Limiting clusters based on acceptable case and control sample sizes and acceptable case/control ratio")
logr::log_print(nclust_size_ratio <- (cl_sizes %>% filter(ratio_josh <=max_ratio,case >= min_sample & control >= min_sample))$cluster %>% print())


# Load in sample files

sample_paths <- paste0(KinshipFlashPCA_path, "flashPCA_lclustering_res_",resolution,"_cluster_",nclust_size_ratio,"_sample.txt")
logr::log_print(sample_paths)
sample_list <- lapply(sample_paths,fread)
# logr::log_print(sample_paths)

# downsample controls to achieve the same case to control ratio as max
downsample_control <- function(sample_list_var){
  case_df <- sample_list_var %>% filter(V6==2)
  ctrl_df_orig <- sample_list_var %>% filter(V6==1)
  logr::log_print(case_df_num <- nrow(case_df))
  logr::log_print(paste0("Original Ctrl number is ", ctrl_df_orig_num <- nrow(ctrl_df_orig)))
  logr::log_print(ctrl_df_desired_num <- round(case_df_num/max_ratio))
  logr::log_print(paste0("Ctrls now downsampled to ", nrow(ctrl_df_downsample <- ctrl_df_orig[sample(nrow(ctrl_df_orig),ctrl_df_desired_num),])))
  return(rbind(case_df, ctrl_df_downsample))
}

if(arguments$max_ratio == "FALSE"){
  sample_list_downsampled <- sample_list
} else {
  sample_list_downsampled <- lapply(sample_list,downsample_control)
}

logr::log_print(sample_list_downsampled)

#' 
check_ratios <- function(list_df){
  logr::log_print(head(list_df))
  temp <- list_df %>% group_by(V6) %>% summarize(count_var = n()) 
  logr::log_print(temp$count_var[2]/temp$count_var[1])
}
lapply(sample_list_downsampled,check_ratios)


sample_downsampled_df <- do.call("rbind",sample_list_downsampled)
cluster_label <- nrow(cl_sizes)



tryCatch({
  logr::log_print(save_path <- paste0(KinshipFlashPCA_path,"flashPCA_lclustering_res_",resolution,"_cluster_",cluster_label,"_min_sample_",min_sample,"_max_ratio_",arguments$max_ratio,"_sample.txt"))

  write.table(x = sample_downsampled_df, file = save_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}, error = function(e) {
  # Handle the error here
  print(paste("An error occurred while saving the new ped file:", e$message))
})

# close log file
logr::log_print("Finished")
logr::log_close()