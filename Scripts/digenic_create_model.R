# digenic_create_model.R: Creates a digenic model out of a single model

'Usage: 
  digenic_create_model.R --resolution_var=<resolution_var> --model=<model>
  
  Options:
  -h --help
  --resolution_var=<resolution_var> resolution for clustering. Typically 0.1-0.4
  --model=<model>
' -> doc
library(tidyverse)
library(data.table)
library(here)
source(here("Scripts","fp_forest_plot_functions.R"))
library(docopt)
arguments <- docopt(doc, version = 'digenic_create_model.R 1.1')

# debug
# arguments <- list()
# arguments$resolution_var <- "0_2"
# arguments$model <- "URFUNC_min80"

#########  ############
dir.create(log_folder <- here("Data","digenic_create_model"))
time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_")
logr::log_open(paste0(log_folder,"/", time_case_prefix,"digenic_create_model.log"))
logr::log_print(arguments)

#########  ############

# load_model <- function(arguments){
#   resolution_var <- arguments$resolution_var
#   model <- arguments$model
#   return(df_list <- cluster_sample_and_genotype_list_fp(paste0(here(),"/"),resolution_var, 0, model, matrix_geno = "matrix"))
# }

create_digenic_model <- function(arguments,atav_df){
  atav_distinct_df <- distinct(atav_df, `Sample Name`,`Gene Name`,.keep_all = TRUE)
  sample_names <- sort(unique(atav_distinct_df$`Sample Name`))
  y <- function(z){
    temp_list <- list()
    for(i in 1:(nrow(z)-1)){
      for(j in (i+1):nrow(z)){
        temp_list[[length(temp_list)+1]]<-data.frame(`Sample Name` = z$`Sample Name`[i],Digene = paste0(z$`Gene Name`[i],"_",z$`Gene Name`[j]), gene1 = z$`Gene Name`[i], gene2 = z$`Gene Name`[j], cluster = z$cluster[1], check.names = TRUE)
      }
    }
    return(temp_list)
  }
  
  # creates a distinct sample/gene data frame for each sample
  distinct_df_list <- lapply(sample_names, function(x) atav_distinct_df %>% filter(`Sample Name` == x) %>% arrange(`Gene Name`))
  # for those sample that have at least 2 QVs in "model", create a new df with columns Sample.Name, 
  digenic_rows_df_list <- lapply(distinct_df_list, function(x) if(nrow(x) > 1) y(x))
  # removes null list
  digenic_rows_df_list_no_null <- digenic_rows_df_list %>% discard(is.null)
  # creates a list of DFs instead of a list of lists
  digenic_rows_df_list_no_null_dfs <- lapply(digenic_rows_df_list_no_null, function(x) do.call(rbind, x))
  # creates a large df
  return(digenic_rows_df <- do.call( rbind, digenic_rows_df_list_no_null_dfs))
}

save_digenic_model <- function(digenic_df_for_collapsing,digenic_df,arguments){
  for(i in unique(digenic_df$cluster)){
    dir.create(save_path <- here("Results","Collapsing",paste0("LClust_res_",arguments$resolution_var,"_cluster_",i,"_FlashColl_07/",arguments$model,"_digenic/")))
    print(save_name <- paste0(save_path,arguments$model,"_digenic_genotypes.csv"))
    print(nrow(save_df <- digenic_df_for_collapsing %>% filter(cluster == i)))
    write_csv( x = save_df,  file =  save_name, col_names = TRUE)
    
    print(save_name <- paste0(save_path,arguments$model,"_digenic_details.csv"))
    write_csv( x = digenic_df %>% filter(cluster == i),  file =  save_name, col_names = TRUE)
  }
}

load_model_cluster <- function(arguments, cluster){
  load_path <- list.files(path = here("Results","Collapsing",paste0("LClust_res_",arguments$resolution_var,"_cluster_",cluster,"_FlashColl_07/",arguments$model)), pattern = "genotypes.csv$", full.names = TRUE)
  atav_df <- fread(load_path)
  atav_df$cluster <- cluster
  return(atav_df)
}

logr::log_print(cl_path <- list.files(path = here("Results","Plots"), pattern = paste0("res_", arguments$resolution_var,"_cluster_sizes.txt$"), full.names = TRUE))
cl_sizes <- fread(cl_path)
nclust <- cl_sizes$cluster
start_list <- lapply(nclust, function(x) load_model_cluster(arguments, x))
start_df <- do.call( rbind, start_list)
nrow(digenic_df <- create_digenic_model(arguments,start_df))

# nrow(distinct(start_df, `Sample Name`,`Gene Name`,.keep_all = TRUE))

digenic_df$`Gene Name` <- digenic_df$gene1
digenic_df$`Sample Name` <- digenic_df$Sample.Name
nrow(start_merge_df <- distinct(start_df, `Sample Name`, `Gene Name`, .keep_all = TRUE))
nrow(end_merge_df <- distinct(merge(digenic_df,start_merge_df, all = TRUE), `Sample Name`, Digene, .keep_all = TRUE) %>% drop_na(Digene))

save_digenic_model(end_merge_df, digenic_df,arguments)

logr::log_close()

