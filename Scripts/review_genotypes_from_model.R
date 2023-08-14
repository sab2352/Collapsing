# review_genotypes_from_model.R: combines genotype files from a model across clusters. saves total, case only, and ctrl only

'Usage: 
  review_genotypes_from_model.R --resolution=<resolution> [--min_sample=<min_sample>] [--case_group=<case_group>] [--combine_ctrl_output=<combine_ctrl_output>] [--debug=<debug>]


  Options:
  -h --help
    --resolution=<resolution> resolution used to create clusters
    --min_sample=<min_sample> [default: NA]
    --case_group=<case_group> [default: temp]
    --combine_ctrl_output=<combine_ctrl_output> [default: False]
    --debug=<debug> [default: FALSE]
' -> doc

library(docopt)
library(tidyverse)
library(data.table)
library(parallel)
library(here)
library(yaml)
"%!in%" <- Negate("%in%")

arguments <- docopt(doc, version = 'review_genotypes_from_model.R')

if(arguments$debug == "TRUE"){
  arguments <- list()
  arguments$min_sample <- "5"
  arguments$resolution <- "0_04"
  arguments$case_group <- "temp"
  arguments$combine_ctrl_output <- "False"
}

#' Initialize log file
#' @param arguments The arguments list with the scripts arguments
initialize_logfile<- function(arguments){
  dir.create(here("Data"))
  dir.create(log_folder <- here("Data","review_genotypes_from_model_Log"))
  time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_")
  logr::log_open(paste0(log_folder,"/", time_case_prefix,"review_genotypes_from_model_logfile.log"))
  logr::log_print(arguments)
}

combine_genotypes <- function(model) {
  files_list <- list.files(path = here("Results","Collapsing",paste0("LClust_res_",arguments$resolution, "_cluster_0_FlashColl_07"),model), pattern = paste0(model,"_genotypes.csv$"),full.names = TRUE)
  #write ifelse if files_list is_empty, use comphet

  if (!is_empty(files_list)) {
    name <- "_genotypes"
  } else {
    name <- "_comphet"
  }
  genotype_name_list <- lapply(analyzed_clusters, 
                              function(x) list.files(path = here("Results","Collapsing",paste0("LClust_res_",arguments$resolution, "_cluster_",x,"_FlashColl_07"),model), pattern = paste0(model, name,".csv$"),full.names = TRUE))
  genotype_list <- mclapply(1:length(genotype_name_list), function(x) fread(genotype_name_list[[x]], header = TRUE, sep = ",", fill = TRUE), mc.cores = 1)
  genotype_df <- do.call(rbind, genotype_list)

  logr::log_print(sprintf("Writing CSV files for %s", model))

  if (arguments$combine_ctrl_output == "True") {
    write.csv(x = genotype_df, file = here("Results","genotype_review",paste0(arguments$case_group, "_model_", model,"_resolution_",arguments$resolution,"_min_sample_",arguments$min_sample ,".csv")), row.names = FALSE)
  }

  if (name == "_genotypes") {
    write.csv(x = genotype_df %>% filter(`Sample Phenotype` =="case"), file = here("Results","genotype_review",paste0(arguments$case_group, "_model_", model,"_resolution_",arguments$resolution,"_min_sample_",arguments$min_sample,"_case_only.csv")), row.names = FALSE)
    if (arguments$combine_ctrl_output == "True") {
      write.csv(x = genotype_df %>% filter(`Sample Phenotype` =="ctrl"), file = here("Results","genotype_review",paste0(arguments$case_group, "_model_", model,"_resolution_",arguments$resolution,"_min_sample_",arguments$min_sample,"_ctrl_only.csv")), row.names = FALSE)
    }
  } else {
    write.csv(x = genotype_df %>% filter(`Sample Phenotype (#1)` =="case"), file = here("Results","genotype_review",paste0(arguments$case_group, "_model_", model,"_resolution_",arguments$resolution,"_min_sample_",arguments$min_sample,"_case_only.csv")), row.names = FALSE)
    if (arguments$combine_ctrl_output == "True") {
      write.csv(x = genotype_df %>% filter(`Sample Phenotype (#1)` =="ctrl"), file = here("Results","genotype_review",paste0(arguments$case_group, "_model_", model,"_resolution_",arguments$resolution,"_min_sample_",arguments$min_sample,"_ctrl_only.csv")), row.names = FALSE)
    }
  }
}

###Main#####
initialize_logfile(arguments)
logr::log_print("log file initialized")

dir.create(save_path <- here("Results","genotype_review"))

yaml_file <- read_yaml(here("Input", "collapsing.yaml"))
model_str <- yaml_file$PERMUTATION_VARIABLE$perm_models
models <- as.list(strsplit(model_str, ",")[[1]])

logr::log_print(sprintf("Cluster size text file is %s",cl_size_path <- list.files(here("Results","Plots"), pattern = paste0("*lclustering_res_", arguments$resolution, "_cluster_sizes.txt$"), full.names = TRUE)))
cl_sizes <- fread(cl_size_path)
ifelse(arguments$min_sample == "NA",analyzed_clusters <- cl_sizes$cluster,analyzed_clusters <- filter(cl_sizes,case >= as.numeric(arguments$min_sample) & control >= as.numeric(arguments$min_sample)  )$cluster)

# # "Folders of included clusters"
# nclust_folders <- lapply(analyzed_clusters, function(x) here("Results","Collapsing",paste0("LClust_res_",arguments$resolution, "_cluster_",x,"_FlashColl_07")))
# # "Models in folders"
# models.list <- lapply(nclust_folders, function(cl_folder) list.dirs(cl_folder, recursive = FALSE, full.names = FALSE))
# models <- Reduce(intersect, models.list)
logr::log_print("models")
logr::log_print(models)

lapply(models, combine_genotypes)
