'Usage: 
  cohortSelection_summary.R [--sample_path_w_filename=<sample_path_w_filename>] [--debug=<debug>] [--case_group=<case_group>] [--exclude_list=<exclude_list>]
  
  Options:
  -h --help
    --sample_path_w_filename=<sample_path_w_filename> path to csv file downloaded from https://sequence.igm.cumc.columbia.edu/ [default: /nfs/goldstein/software/atav_home/data/sample/igm_lims_sample_master_list.tsv]
    --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
    --exclude_list=<exclude_list> [default: NA]
    --debug=<debug> [default: FALSE]
' -> doc

library(data.table)
library(here)
library(tidyverse)
library(logr)
library(docopt)
"%!in%" <- Negate("%in%")
arguments <- docopt(doc, version = "regenieCovPhenoSplit.R 1.1")

if(arguments$debug == "TRUE"){
  arguments <- list()
  arguments$case_group <- "debug"
  arguments$sample_path_w_filename <- here("Data","sample_result_16940026_1687870417.csv")
  arguments$exclude_list <- c("IGM-GHARMaGICP22uu18884625,IGM-GHARMaGIC22uuSF086-IV-06-P2fSF086,IGM-GHARMaGICP22uuA111fSF178,IGM-GHARMaGICP22uuA048fSF343,IGM-GHARMaGIC22uuSF139-III-03-P2fSF139,IGM-GHARMaGIC22uuSF184-II-01-P2fSF184,IGM-GHARMaGIC22uuSF036-II-1-P2fSF036,IGM-GHARMaGIC22uuSF310-II-02-P2fSF310")
                              
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
  initialize_logfile(arguments, "cohortSelection_summary")
}, error = function(e) {
  # Handle the error here
  print(paste("An error occurred while creating the logfile:", e$message))
})

tryCatch({
  case_included_df <- fread(case_included <- list.files(path = here("Data"), pattern = "*case.seq.csv$", full.names = TRUE))
  case_excluded_df <- fread(case_excluded_df <- list.files(path = here("Data"), pattern = "\\d_excluded_cases.csv$", full.names = TRUE))
  seq_df <- fread(arguments$sample_path_w_filename)
}, error = function(e) {
  # Handle the error here
  print(paste("An error occurred while while loading in csv files:", e$message))
})
tryCatch({
  analyzed_families <- unique(c(case_excluded_df$FamilyID, case_included_df$FamilyID))
  logr::log_print(sprintf("There are %i analyzed families", length(analyzed_families)))
  seq_df_incl_families <- seq_df %>% filter(FamilyID %in% analyzed_families, sample_status == "In DragenDB", sample_internal_name %!in% unlist(strsplit(arguments$exclude_list,",")))
  family_structure <- seq_df_incl_families %>% group_by(FamilyID) %>% summarize(family_size = n())
  family_structure_cts <- family_structure %>% group_by(family_size) %>% summarize(count_var = n())
  logr::log_print("Family sizes are as follows")
  logr::log_print(family_structure_cts)
  mult_case_fms <- rbind(case_included_df, case_excluded_df) %>% filter(sample_internal_name %!in% unlist(strsplit(arguments$exclude_list,","))) %>% group_by(FamilyID) %>% summarize(family_size = n()) %>% group_by(family_size) %>% summarize(count_var = n()) 
}, error = function(e) {
  # Handle the error here
  print(paste("An error occurred while summarizing the family structure:", e$message))
})

tryCatch({
  logr::log_print("Finished")
  logr::log_close()
}, error = function(e) {
  # Handle the error here
  print(paste("An error occurred while closing log file:", e$message))
})
