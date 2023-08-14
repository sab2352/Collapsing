# table_sample_counts.R: creates a table with cases and controls at each step

'Usage: 
  table_sample_counts.R --resolution_var=<resolution_var> --min_sample=<min_sample> [--case_group=<case_group>]
  
  Options:
  -h --help
  --min_sample=<min_sample> minimum number of cases/controls of included clusters
  --resolution_var=<resolution_var> resolution for clustering.
  --debug=<debug> [default: FALSE]
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]

' -> doc

library(docopt)
library(tidyverse)
library(data.table)
library(here)
library(logr)
library(yaml)
library(sjPlot)
"%!in%" <- Negate("%in%")
arguments <- docopt(doc, version = 'table_sample_counts.R 1.1')

if(arguments$debug == "TRUE"){
  arguments <- list()
  arguments$min_sample <- "15"
  arguments$resolution_var <- "0_5"
  arguments$case_group <- "debug"
}

print_pre_post <- function(name_of_step, path_var){
    
  if(length(path_var) > 0 && file.exists(path_var)){
    ped_df_var <- fread(path_var)
    logr::log_print(sprintf("%s ped contains %i cases and %i controls", name_of_step,
                            case_num <- nrow(ped_df_var %>% filter(V6==2)),
                            ctrl_num <- nrow(ped_df_var %>% filter( V6==1))))
    return(c(case_num, ctrl_num))
  } else {
    logr::log_print(sprintf("%s ped does not exist", name_of_step))
    return(c(NA, NA))
  }

}

print_steps_from_root <-function(arguments){
  resolution_var <- arguments$resolution_var
  min_sample <- as.numeric(arguments$min_sample)
  yaml_file <- read_yaml(here("Input", "collapsing.yaml"))
  
  # original case list
  case_list <- fread(here("Input",yaml_file$USER_VARIABLE$case_path), header = FALSE)
  logr::log_print(sprintf("Originally case list had %i samples", nrow(case_list)))
  case_list_matrix <- data.frame("Analysis Step" = rep("0 Case List",2),
                              "Case-Ctrl" = c("Case","Control"),
                              "Count" = c(nrow(case_list), 0))
  
  # **********
  # Cohort-----
  # **********
  excl_case_seq_path <- list.files(here("Data"), pattern = "\\d_excluded_cases.csv$", full.names = TRUE)
  if(file.exists(excl_case_seq_path)){
    excl_case_seq <- fread(excl_case_seq_path)
  } else {
    logr::log_print("There are no excluded cases")
    return(c(NA, NA))
  }

  # included and excluded
  case_seq_path <- list.files(here("Data"), pattern = "*case.seq.csv$", full.names = TRUE)
  case_seq <- fread(case_seq_path)
  logr::log_print(sprintf("Originally number of samples is %i of which %i are included and %i are excluded", nrow(case_seq) + nrow(excl_case_seq) ,nrow(case_seq), nrow(excl_case_seq)))

  print_pre_post_loop_fxn<-function(x,input_list){
    ped_path <- list.files(input_list[1], pattern = input_list[2], full.names = TRUE)
    case_ctrl_nums <- print_pre_post(input_list[3], ped_path)
    return_matrix <- data.frame("Analysis Step" = rep(paste0(x, " ",input_list[3]),2),
                                "Case-Ctrl" = c("Case","Control"),
                                "Count" = case_ctrl_nums)
    return(return_matrix)
  }
  
  input_list_var <- list(c(here("Data"), "*ped.txt$", "Pre Sample Check"),
       c(here("Results","SampleFileCheck"), "*existing.sample.txt$", "Post Sample Check"),
       c(here("Results","Coverage"), "*existing.sample.txt$","Initial Coverage"),
       c(here("Results","KinshipFlashPCA"), "*_KinshipFlashPCA_existing.sample.txt$", "PCA Generation"),
       c(here("Results","KinshipFlashPCA"), "*KinshipFlashPCA_kinship_pruned_sample.txt$", "Kinship Pruned"),
       c(here("Results","KinshipFlashPCA"), paste0("*flashPCA_lclustering_res_", yaml_file$COLLAPSING_VARIABLE$resolution,
                                "_cluster_",yaml_file$COLLAPSING_VARIABLE$mixed_ancestry_cluster,"_min_sample_",yaml_file$COLLAPSING_VARIABLE$min_sample,"_max_ratio_",
                                gsub("\\.","_", yaml_file$COLLAPSING_VARIABLE$max_ratio),"_sample.txt"), "Mixed Ancestry Initial"),
       c(paste0(here("Results","Coverage"), "/CoverageFlashLClust_res_", yaml_file$COLLAPSING_VARIABLE$resolution,
                "_cluster_",yaml_file$COLLAPSING_VARIABLE$mixed_ancestry_cluster,"_min_sample_",yaml_file$COLLAPSING_VARIABLE$min_sample,"_max_ratio_",
                gsub("\\.","_", yaml_file$COLLAPSING_VARIABLE$max_ratio),"_07_mixedCluster"), "*_existing.sample.txt$","Mixed Ancestry Coverage" ),
       c(paste0(here("Results","Collapsing"), "/LClust_res_", yaml_file$COLLAPSING_VARIABLE$resolution,
                "_cluster_",yaml_file$COLLAPSING_VARIABLE$mixed_ancestry_cluster,"_min_sample_",yaml_file$COLLAPSING_VARIABLE$min_sample,"_max_ratio_",
                gsub("\\.","_", yaml_file$COLLAPSING_VARIABLE$max_ratio),"_mixedCluster_FlashColl_07/dominantNoneMAF"),"*_existing.sample.txt$","Mixed Ancestry Variant Calling"))
  
  df_list <- lapply(1:length(input_list_var),function(x) print_pre_post_loop_fxn(x,input_list_var[[x]]))
  
  logr::log_print(cohort_matrix_pre <- do.call(rbind, df_list))
  cohort_matrix <- rbind(case_list_matrix, cohort_matrix_pre)
  # **********
  # Write cohort and mixed ancestry-----
  # **********
  logr::log_print(write_df <- cohort_matrix %>% spread(Case.Ctrl, Count))
  tab_df(write_df, file = save_path <- here(output_folder, paste0(time_case_prefix,"cohort_counts.doc")), alternate.rows=TRUE, col.header = names(write_df),use.viewer = TRUE,
         show.footnote = TRUE, footnote = arguments$footnote)
  write_csv(x = write_df, file = save_path <- here(output_folder, paste0(time_case_prefix,"cohort_counts.csv")))#, col_names = TRUE, quote = FALSE)
  logr::log_print(sprintf("Finished writing cohort table to %s",save_path ))
  
  
  # **********
  # Clusters-----
  # **********
  # # load in original cluster sizes
  logr::log_print(sprintf("Cluster size text file is %s",cl_size_path <- list.files(here("Results","Plots"), pattern = paste0("*lclustering_res_", resolution_var, "_cluster_sizes.txt$"), full.names = TRUE)))
  cl_sizes <- fread(cl_size_path)
  cluster_matrix <- data.frame(Clusters = cl_sizes$cluster)
  logr::log_print(sprintf("Original clusters were the following"))
  logr::log_print(cl_sizes)
  logr::log_print(sprintf("Analyzed clusters contain %i cases and %i controls", 
                          sum((cl_sizes %>% filter(case >= min_sample & control >= min_sample))$case),
                          sum((cl_sizes %>% filter(case >= min_sample & control >= min_sample))$control)  ))
  cluster_matrix$Analyzed <- "No"
  analyzed_clusters <- (cl_sizes %>% filter(case >= min_sample & control >= min_sample))$cluster
  cluster_matrix$Analyzed[analyzed_clusters+1] <- "Yes"
  cluster_matrix$`Initial Case` <- cl_sizes$case
  cluster_matrix$`Initial Ctrl` <- cl_sizes$control
  # was PCA stuff run
  y<-function(x){
    list.files(here(paste0("Results/KinshipFlashPCA/PCAFlashLClust_res_",resolution_var,"_cluster_",x,"_07")), pattern = "*lashpca_pruned_sample_file.txt$", full.names = TRUE)
  }
  sample_file_location <- lapply(cl_sizes$cluster,y)
  summarize_sample_file <- function(x,sample_file_paths){
    if(length(sample_file_paths[[x]])>0){
      sample_file_var <- fread(sample_file_paths[[x]])
      case_num <- nrow(sample_file_var %>% filter(V6==2))
      ctrl_num <- nrow(sample_file_var %>% filter(V6==1))
      return(c(case_num, ctrl_num))
    } else {
      logr::log_print(paste0("No file for cluster",x-1))
      return(c(NA, NA))
    }
  }
  case_ctrl_nums <- lapply(cl_sizes$cluster+1,summarize_sample_file,sample_file_location)
  cluster_matrix$`PCA Ancestry Excl Case` <- sapply(cl_sizes$cluster+1, function(x) case_ctrl_nums[[x]][1])
  cluster_matrix$`PCA Ancestry Excl Ctrl` <- sapply(cl_sizes$cluster+1, function(x) case_ctrl_nums[[x]][2])
  
  # was master genotype file generated
  y<-function(x){
    list.files(here(paste0("Results/Collapsing/LClust_res_",resolution_var,"_cluster_",x,"_FlashColl_07/dominantNoneMAF")), pattern = "*_existing.sample.txt$", full.names = TRUE)
  }
  sample_file_location <- lapply(cl_sizes$cluster,y)
  case_ctrl_nums <- lapply(cl_sizes$cluster+1,summarize_sample_file,sample_file_location)
  cluster_matrix$`Variant Call Case` <- sapply(cl_sizes$cluster+1, function(x) case_ctrl_nums[[x]][1])
  cluster_matrix$`Variant Call Ctrl` <- sapply(cl_sizes$cluster+1, function(x) case_ctrl_nums[[x]][2])
  logr::log_print(cluster_matrix)
  write_csv(x = cluster_matrix, file = save_path <- paste0(output_folder, "/",time_case_prefix,"cluster_counts.csv"))#, col_names = TRUE, quote = FALSE)
  logr::log_print(sprintf("Finished writing cluster table to %s",save_path ))
  

}

#' Initialize log file
#' @param arguments The arguments list with the scripts arguments
initialize_logfile <- function(arguments, label_var) {
  dir.create(here("Data"))
  dir.create(log_folder <- here("Data", paste0(label_var,"_Log")))
  logr::log_open(here(log_folder, paste0(time_case_prefix, label_var,"_logfile.log")))
  logr::log_print(arguments)
}

# **********
# Main-----
# **********
dir.create(here("Data"))
dir.create(here("Results"))
dir.create(output_folder <- paste0(here("Results","SampleTracking")))
time_case_prefix <- paste0(gsub(":", "_", gsub(" ", "_", gsub("-", "_", Sys.time()))), "_", arguments$case_group, "_")

# Create log files
tryCatch({
  initialize_logfile(arguments, "SampleTracking")
}, error = function(e) {
  # Handle the error here
  print(paste("An error occurred while creating the logfile:", e$message))
})

tryCatch({
  print_steps_from_root(arguments)
}, error = function(e) {
  # Handle the error here
  print(paste("An error occurred in everything else:", e$message))
})

tryCatch({
  logr::log_print("Finished")
  logr::log_close()
}, error = function(e) {
  # Handle the error here
  print(paste("An error occurred while closing log file:", e$message))
})
