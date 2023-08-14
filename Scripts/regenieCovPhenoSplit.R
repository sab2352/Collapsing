'Usage: 
  regenieCovPhenoSplit.R --regenie_folder=<regenie_folder> --pcs=<pcs> [--case_group=<case_group>] [--debug=<debug>]
  
  Options:
  -h --help
    --regenie_folder=<regenie_folder> folder with regenie project
    --pcs=<pcs> is the allele frequencies used for masking (includes singleton as default)
    --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
    --debug=<debug> [default: FALSE]
  #eventually will want the ability to do multiple phenotypes in the pheno file
  #best would be to use sample names to create phenotypes here - probably utilizing a new folder in /Regenie/
  
' -> doc

library(data.table)
library(here)
library(dplyr)
library(logr)
library(docopt)
arguments <- docopt(doc, version = "regenieCovPhenoSplit.R 1.1")

if(arguments$debug == "TRUE"){
  arguments <- list()
  arguments$regenie_folder <- here("Results","Regenie","0_4_cluster_11_min_sample_100_max_ratio_True")
  arguments$case_group <- "debug"
  arguments$pcs <- "10"
}

time_case_prefix <- paste0(gsub(":", "_", gsub(" ", "_", gsub("-", "_", Sys.time()))), "_", arguments$case_group, "_")
#' Initialize log file
#' @param arguments The arguments list with the scripts arguments
initialize_logfile <- function(arguments, label_var) {
  dir.create(here("Data"))
  dir.create(log_folder <- here("Data", paste0(label_var,"_Log")))
  logr::log_open(here(log_folder, paste0(time_case_prefix, label_var,"_logfile.log")))
  logr::log_print(arguments)
}

# Create log files
tryCatch({
  initialize_logfile(arguments, "regenieCovPhenoSplit")
}, error = function(e) {
  # Handle the error here
  print(paste("An error occurred while creating the logfile:", e$message))
})

# rest of code
tryCatch({
  phenotypes<-fread(here(arguments$regenie_folder,paste0("CovPheno_pc", arguments$pcs),"phenotypes.txt"))

  # for every column in the "phenotypes" dataframe (or list), it determines the class (i.e., the data type) of the column and stores these classes in the "colclass" vector. For example, if "phenotypes" was a dataframe with three columns: "height" of class numeric, "weight" of class numeric, and "hair_color" of class character, sapply(phenotypes,class) would return a vector: ("numeric", "numeric", "character"). This vector is then stored in "colclass". This kind of operation might be useful if you wanted to quickly get an overview of the data types in a dataframe.
  sapply(phenotypes,class)->colclass
  
  phenotypesBT<-cbind(phenotypes[,1:2],select(phenotypes,which(colclass=="integer")))
  
  phenotypesQT<-cbind(phenotypes[,1:2],select(phenotypes,which(colclass=="numeric")))
  
  write.table(phenotypesBT,here(arguments$regenie_folder,paste0("CovPheno_pc", arguments$pcs),"phenotypesBT.txt"),row.names = F,sep="\t",quote = F)
  write.table(phenotypesQT,here(arguments$regenie_folder,paste0("CovPheno_pc", arguments$pcs),"phenotypesQT.txt"),row.names = F,sep="\t",quote = F)
}, error = function(e) {
  # Handle the error here
  print(paste("An error in rest of code:", e$message))
})

# close log file
logr::log_print("Finished")
logr::log_close()