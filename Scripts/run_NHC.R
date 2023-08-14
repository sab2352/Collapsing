# run_NHC.R: runs NHC (https://github.com/casanova-lab/NHC)

'Usage: 
  run_NHC.R --cluster=<cluster> --resolution_var=<resolution_var> --model=<model> [--path_to_nhc_git=<path_to_nhc_git>] [--gene_list_case_path=<gene_list_case_path>] [--gene_list_ctrl_path=<gene_list_ctrl_path>]  [--min_sample=<min_sample>] [--case_group=<case_group>] [--max_ratio=<max_ratio>] [--num_pcas=<num_pcas>] [--mixed_ancestry_logical=<mixed_ancestry_logical>] [--w=<w>] [--b=<b>] [--m=<m>]
  
  Options:
  -h --help
  --cluster=<cluster> define cluster to run NHC on
  --path_to_nhc_git=<path_to_nhc_git> define path to NHC git hub repo. See https://github.com/casanova-lab/NHC for NHC documentation [default: /home/jm4279/old_goldstein_lab_folder/jm4279/github_repos/NHC/]
  --model=<model> will run on a model if specified. If specified, no need for gene_list_case_path or gene_list_ctrl_path
  --gene_list_case_path=<--gene_list_case_path> path to list of cases and genes harboring QV. Best practice is for this list to live in a folder /Results/NHC [default: NA]
  --resolution_var=<resolution_var> resolution for clustering. Typically 0.1-0.4
  --gene_list_ctrl_path=<--gene_list_ctrl_path> path to list of controls and genes harboring QV. [default: NA]
  --case_group=<case_group> [default: temp_group]
  --min_sample=<min_sample> minimum size of case/control allowed to include cluster into combined cluster. Need to overide default if using combined cluster. Otherwise, will  [default: NA]
  --max_ratio=<max_ratio> this is added to allow loading of combined ancestry cluster [default: NA]
  --mixed_ancestry_logical=<mixed_ancestry_logical> [default: FALSE]
  --b=<b> Hub gene removal is to avoid giant clusters that are formed due to the large number of interactions with hub genes. The connectivity of each gene is determined by the number of PPIs above STRING score 0.9 (Data_Network_Connectivity.txt). The default value (-b 100) means: skipping the genes having more than 100 PPIs with edge-weight>0.9 for clustering. If users want to include all genes for clustering, use (-b 0). [default: 100]
  --m=<m> merge overlapped clusters (overlapping ratio = common/union genes) The initial outputted gene clusters are then iteratively merged, if one is a superset/subset of another, or the two most-overlapping clusters sharing over 50% (default) genes, thereby reducing the number of gene clusters that are more distinct from each other. [default: 0.5]
  --w=<w> float	edge-weight cutoff, based on STRING score 0.7~1	0.99 Stringent edge-weight cutoff is used to converge the gene clusters of the highest biological relevance. If the case cohort is small or the gene candidates are few, then users could relax the edge-weight cutoff to 0.95 or 0.9, but no lower than 0.7 (as STRING determines 0.7 as confidence cutoff). [default: 0.99]
  --debug=<debug> [default: FALSE]

' -> doc
library(tidyverse)
library(data.table)
library(here)
library(docopt)
library(logr)
arguments <- docopt(doc, version = 'run_NHC.R 2.2')
# debugging
if(arguments$debug == "TRUE"){
  arguments <- list()
  arguments$cluster <- "11"
  arguments$resolution_var <- "0_4"
  arguments$model <- "URFUNC_min80_igm"
  arguments$min_sample <- "40"
  arguments$max_ratio <- "True"
  arguments$mixed_ancestry_logical <- "True"
  arguments$gene_list_case_path <- "NA"
  # arguments$path_to_nhc_git <- "W:\\old_goldstein_lab_folder\\jm4279\\github_repos\\NHC\\"
  arguments$path_to_nhc_git <- "/mnt/w/old_goldstein_lab_folder/jm4279/github_repos/NHC/"
  arguments$w <- "0.99"
  arguments$m <- "0.5"
  arguments$b <- "100"
} 

# Create Logfile----
tryCatch({
  dir.create(log_folder <- here( "Data"))
  dir.create(log_folder <- here( "Data","run_NHC_log"))
  time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_", arguments$model,"_","mixedancestry_",arguments$mixed_ancestry_logical,"_")
  logr::log_open(log_folder,paste0(time_case_prefix,"run_NHC_logfile.log"))
  logr::log_print(arguments)
  logr::log_print(R.version)
}, error = function(e) {
  message("Caught an error making log: ", e$message)
})
# create directories
NHC_path <- here("Results","NHC")
if(!dir.exists(NHC_path)) {
  dir.create(NHC_path)
} else {
  logr::log_print(paste0("Directory ", NHC_path, " already exists"))
}
output_directory <- here("Results","NHC", paste0(time_case_prefix,"output"))
if(!dir.exists(output_directory)) {
  dir.create(output_directory)
} else {
  logr::log_print(paste0("Directory ", output_directory, " already exists"))
}


# Load in PCAs----
tryCatch({
  if(arguments$mixed_ancestry_logical == "FALSE"){
    pc_path <- here("Results","KinshipFlashPCA",paste0("PCAFlashLClust_res_", arguments$resolution_var,"_cluster_",arguments$cluster, "_07"))
  } else {
    pc_path <- here("Results","KinshipFlashPCA",paste0("PCAFlashLClust_res_", arguments$resolution_var,"_cluster_",arguments$cluster, "_min_sample_",arguments$min_sample,"_max_ratio_",arguments$max_ratio, "_07_mixedCluster"))
  }
    
  logr::log_print(sprintf("This is the pc_path %s", pc_path))#used to have dir here but removed for now
  logr::log_print(evf <- list.files(pc_path, "*eigenvectors.txt$", full.names = TRUE))
  # logr::log_print(evf <- list.files(pc_path, "*_flashpca_eigenvectors$", full.names = TRUE))
  logr::log_print(sprintf("This is the evf (path to the PCs) %s", evf))
  # reformat PCAs
  logr::log_print("Now reading in PCs")
  pcs <- fread(here(evf))
  logr::log_print("Finished reading in PCs")
}, error = function(e) {
  message("Caught an error loading in PCs: ", e$message)
})

tryCatch({
  pcs_formatted <- pcs %>% dplyr::select(sampleID = FID, pca_1 = U1, pca_2 = U2, pca_3 = U3)
  logr::log_print("Now writing reformatted in PCs")
  write.table(x = pcs_formatted,file = pcs_file <- here(output_directory,paste0(time_case_prefix, "PCAs.txt")),quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE )
  logr::log_print("Finished writing reformatted in PCs")
}, error = function(e) {
  message("Caught an error reformatting and writing: ", e$message)
})


# load in case genes
if(arguments$gene_list_case_path != "NA"){
  logr::log_print("arguments$gene_list_case_path != NA so reading in case_path %s",arguments$gene_list_case_path)
  gene_list_case <- fread(arguments$gene_list_case_path, header = TRUE) #
  logr::log_print("Finished reading %s",arguments$gene_list_case_path)
} else {
  logr::log_print("arguments$gene_list_case_path == NA")
  if(arguments$mixed_ancestry_logical == "FALSE"){
    logr::log_print("Not ancestry cluster")
    collapsing_folder <- here("Results","Collapsing",paste0("LClust_res_",arguments$resolution_var,"_cluster_",arguments$cluster, "_FlashColl_07/",arguments$model))
  } else {
    logr::log_print("Yes ancestry cluster")
    collapsing_folder <- here("Results","Collapsing",paste0("LClust_res_",arguments$resolution_var,"_cluster_",arguments$cluster, "_min_sample_",arguments$min_sample,"_max_ratio_",gsub("\\.","_", arguments$max_ratio), "_mixedCluster_FlashColl_07/",arguments$model))
    logr::log_print(sprintf("Ancestry collapsing folder is %s",collapsing_folder))
  }

  sample_path <- list.files(path = collapsing_folder, "*existing.sample.txt$", full.names = TRUE)[1]
  sample_df <- fread(sample_path)
  logr::log_print("Loaded in Sample File")
  matrix_path <- list.files(collapsing_folder, "*_matrix.txt.gz$", full.names = TRUE)[1]
  matrix_var <- fread(matrix_path)
  logr::log_print("Loaded in Matrix")
  
  matrix_var[,1]
  # 
  x<-function(i,data_matrix){
    sample_name <- names(data_matrix)[i]
    one_zero_sample <- data_matrix[,..sample_name]
    one_zero_gene <- data_matrix[,`sample/gene`]
    genes_positive <- one_zero_gene[one_zero_sample==1]
    # sample_name <- names(data_matrix)[i]
    # genes_positive <- data_matrix[data_matrix[,i]==1,1]
    # if(nrow(genes_positive)>0) return(data.frame("Sample.Name" = rep(sample_name,nrow(genes_positive)), "Gene.Name" = paste0("'",genes_positive$`sample/gene`,"'")))
    if(length(genes_positive)>0) return(data.frame("ID" = sample_name, "GeneList" = paste(as.character(genes_positive), collapse = ",")))
    
    return(NULL)
  }
  logr::log_print("extracting gene sample matrix")
  data_frame_list <- lapply(2:length(names(matrix_var)),x,matrix_var)
  # data_frame_list <- lapply(2,x,matrix_var)
  atav_df <- do.call(rbind, data_frame_list)
  # logr::log_print(atav_df)
  
  case_names <- (sample_df %>% filter(V6 == 2))$V1
  ctrl_names <- (sample_df %>% filter(V6 == 1))$V1
  
  case_df <- atav_df %>% filter(ID %in% case_names)
  logr::log_print(case_df)
  logr::log_print(arguments$gene_list_case_path <- paste0(output_directory, "case_matrix_",time_case_prefix, "case_list.txt"))
  logr::log_print("Writing case matrix")
  write.table(x = case_df, file = arguments$gene_list_case_path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  logr::log_print("Writing ctrl matrix")
  ctrl_df <- atav_df %>% filter(ID %in% ctrl_names)
  write.table(x = ctrl_df, file = arguments$gene_list_ctrl_path <- paste0(output_directory, "ctrl_matrix_",time_case_prefix, "ctrl_list.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

tryCatch({
  logr::log_print(exec_str_case_ctrl <- paste0("python ",arguments$path_to_nhc_git,"NHC_case_control.py -case ",arguments$gene_list_case_path," -ctl ", arguments$gene_list_ctrl_path, " -pc ", pcs_file, " -w ",arguments$w, " -b ",arguments$b, " -m ",arguments$m,  " -o ", here(output_directory,paste0(time_case_prefix,"case_ctrl_output.txt"))))
  
  setwd(arguments$path_to_nhc_git)
  logr::log_print("Running case control")
  system(exec_str_case_ctrl)
  logr::log_print("Finished case control")
}, error = function(e) {
  message("Caught an error running case control: ", e$message)
})

tryCatch({
  logr::log_print(exec_str_case_only <- paste0("python ",arguments$path_to_nhc_git,"NHC_case_only.py -case ",arguments$gene_list_case_pat, " -o ",here(output_directory, paste0(time_case_prefix, "case_only_output.txt")), " -m ",arguments$m, " -w ", arguments$w," -b ",arguments$b)) # 

  setwd(arguments$path_to_nhc_git)
  logr::log_print("Running case only")
  system(exec_str_case_only)
  logr::log_print("Finished case only")
}, error = function(e) {
  message("Caught an error running case only: ", e$message)
})
# Close Logfile----
logr::log_print("Finished NHC Analysis")
logr::log_close()

