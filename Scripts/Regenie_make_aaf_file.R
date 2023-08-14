# Input: AnnoMerged_prefiltered.txt, filtered_genos.csv
# Output: AnnoMerged.txt
'Usage: 
  Regenie_make_aaf_file.R --regenie_folder=<regenie_folder> [--case_group=<case_group>] [--debug=<debug>]
  
  Options:
  -h --help
    --regenie_folder=<regenie_folder> folder with regenie project
    --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
    --debug=<debug> [default: FALSE]
  #eventually will want the ability to do multiple phenotypes in the pheno file
  #best would be to use sample names to create phenotypes here - probably utilizing a new folder in /Regenie/
  
' -> doc

library(docopt)
library(tidyverse)
library(data.table)
library(broom)
library(here)
library(logr)
"%!in%" <- Negate("%in%")
arguments <- docopt(doc, version = 'regenieCovPheno_update.R 1.1')

if(arguments$debug == "TRUE"){
  arguments <- list()
  arguments$regenie_folder <- here("Results","Regenie","0_4_cluster_11_min_sample_100_max_ratio_True")
  arguments$case_group <- "debug"
}


tryCatch({
  ######### Initiate Log ############
  dir.create(log_folder <- here("Data","regenieCovPheno_update"))
  time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_",arguments$case_group, "_")
  logr::log_open(here(log_folder, paste0(time_case_prefix,"regenieCovPheno_update.log")))
  logr::log_print(arguments)
  
}, error = function(e) {
  message("Caught an error making log: ", e$message)
})
tryCatch({
  ####
  logr::log_print("Reading in variants")
  variants <- fread(here(arguments$regenie_folder,"allQCBAFilter","genotypes","filtered_genos.csv")) #,sep = ",",drop=c(10,11,12,21:145)))
  write_df <- data.frame(V1 =variants$`Variant ID`,
                         V2 = variants$`gnomAD Exome global_AF`)
  write_df$V2[is.na(write_df$V2)] <- 0
  logr::log_print("Writing in variants")
  write.table(file = here(arguments$regenie_folder,"allQCBAFilter","genotypes","aaf_file.txt"), x = write_df, quote = FALSE,row.names = FALSE,col.names = FALSE, sep = "\t")
}, error = function(e) {
  message("Caught an error reading variants or writing aaf file: ", e$message)
})

tryCatch({
  logr::log_print("Finished")
  logr::log_close()
}, error = function(e) {
  message("Caught an error closing log: ", e$message)
})
