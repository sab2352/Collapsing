# Input: AnnoMerged_prefiltered.txt, filtered_genos.csv
# Output: AnnoMerged.txt
'Usage: 
  regenieCovPheno_update.R --regenie_folder=<regenie_folder> [--case_group=<case_group>] [--debug=<debug>]
  
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
  anno <- fread(here(arguments$regenie_folder,"Anno","AnnoMerged_prefiltered.txt"),header = F)
  variants <- fread(here(arguments$regenie_folder,"allQCBAFilter","genotypes","filtered_genos.csv")) #,sep = ",",drop=c(10,11,12,21:145)))
  
  MPC3<-variants[`MPC`>=3]
  MPC2<-variants[`MPC`>=2 & `MPC`<3]
  MPC0<-variants[Effect %in% c("missense_variant+splice_region_variant", "missense_variant") & (`Variant ID` %!in% c(MPC2$`Variant ID`, MPC3$`Variant ID`))]
  # missenseMulti<-variants[REVEL>0.35 | `subRVIS Exon Percentile` <50 | PrimateAI >0.50 | `Polyphen Humvar Prediction`!="benign"]
  
  mutate(anno, V3 = case_when(V1 %in% MPC3$`Variant ID`~"MPC3",V1 %in% MPC2$`Variant ID`~"MPC2", V1 %in% MPC0$`Variant ID`~"MPC0", TRUE ~ V3))->anno_updated
  
  # write over the old anno
  write.table(anno_updated,here(arguments$regenie_folder,"Anno","AnnoMerged.txt"),sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
}, error = function(e) {
  message("Caught an error: ", e$message)
})
