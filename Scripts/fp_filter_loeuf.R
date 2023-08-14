# fp_filter_loeuf.R: filters on loeuf

'Usage: 
  fp_filter_loeuf.R --model=<model> --resolution=<resolution> --cluster=<cluster> --loeuf=<loeuf> [--debug=<debug>]

  Options:
  -h --help
    --model=<model> model_to
    --resolution=<resolution> resolution used to create clusters
    --cluster=<cluster> 
    --loeuf=<loeuf>
    --debug=<debug> [default: FALSE]
' -> doc

library(docopt)
library(here)
library(tidyverse)
library(data.table)
"%!in%" <- Negate("%in%")

arguments <- docopt(doc, version = 'fp_filter_loeuf.R 1.1')

if(arguments$debug == "TRUE"){
  arguments <- list()
  arguments$model <- "URPTV_min80"
  arguments$resolution <- "0_2"
  arguments$cluster <- "8"
  arguments$loeuf <- "0.16"
}

#' Initialize log file
#' @param arguments The arguments list with the scripts arguments
initialize_logfile<- function(arguments){
  dir.create(here("Data"))
  dir.create(log_folder <- here("Data","fp_filter_loeuf_Log"))
  time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_")
  logr::log_open(paste0(log_folder,"/", time_case_prefix,"fp_filter_loeuf_logfile.log"))
  logr::log_print(arguments)
}

initialize_logfile(arguments)
logr::log_print("log file initialized")
genotype_name <- list.files(path = here("Results","Collapsing",paste0("LClust_res_",arguments$resolution, "_cluster_",arguments$cluster,"_FlashColl_07"),arguments$model), pattern = paste0(arguments$model,"_genotypes.csv$"),full.names = TRUE)
nrow(genotype_df <- fread(genotype_name))

nrow(genotype_loeuf_df <- genotype_df %>% filter(`gnomAD oe_lof_upper` <= as.numeric(arguments$loeuf)))

dir.create(save_path <- here("Results","Collapsing",paste0("LClust_res_",arguments$resolution, "_cluster_",arguments$cluster,"_FlashColl_07/",arguments$model,"_loeuf_",gsub("\\.","_",sprintf("%.3f",as.numeric(arguments$loeuf))))))


if(nrow(genotype_loeuf_df) == 0){
  logr::log_print("There are 0 QVs in the filtered file")
  fwrite(x = data.frame(), file =  paste0(save_path, "/",arguments$model,"_loeuf_",gsub("\\.","_",sprintf("%.3f",as.numeric(arguments$loeuf))),"_matrix.txt.gz"), compress = "gzip")
}

fwrite(x = genotype_loeuf_df, file = paste0(save_path, "/",arguments$model,"_loeuf_",gsub("\\.","_",sprintf("%.3f",as.numeric(arguments$loeuf))),"_genotypes1.csv"), col.names = TRUE)
logr::log_print("finished")
logr::log_close()
