# model_case_sample.R: lists cases who have variants present in models

'Usage: 
  model_case_sample.R --resolution_var=<resolution_var> --min_sample=<min_sample> [--debug=<debug>]
  
  Options:
  -h --help
  --min_sample=<min_sample> minimum number of cases/controls of included clusters
  --resolution_var=<resolution_var> resolution for clustering.
  --debug=<debug> [default: FALSE]

' -> doc

library(docopt)
library(tidyverse)
library(dplyr)
library(data.table)
library(here)
library(logr)
library(readxl)
library(openxlsx)
arguments <- docopt(doc, version = 'model_case_sample.R 1.1')

if(arguments$debug == "TRUE"){
  arguments <- list()
  arguments$min_sample <- "5"
  arguments$resolution_var <- "0_05"
}

get_cases <- function(model) {
  genes <- fread(here("Results","CMH", paste0(model, "_resolution_", arguments$resolution_var, "_min_sample_", arguments$min_sample, "_top_10_genes.txt")), header = FALSE)
  case_csv <- list.files(here("Results","genotype_review"), pattern= paste0("*", model,"_resolution_",arguments$resolution_var,"_min_sample_",arguments$min_sample,"_case_only.csv"), full.names = TRUE)
  cases <- fread(case_csv)

  if ("Gene Name (#1)" %in% colnames(cases)) {
    genes[["Gene Name (#1)"]] <- paste0("'",genes$V1,"'")
    cases_gene <- cases %>% filter(cases[["Gene Name (#1)"]] %in% genes[["Gene Name (#1)"]])
  } else {
    genes[["Gene Name"]] <- paste0("'",genes$V1,"'")
    cases_gene <- cases %>% filter(cases[["Gene Name"]] %in% genes[["Gene Name"]])
  }

  write.csv(cases_gene, file = here("Results","CMH",model,paste0(model,"_resolution_",arguments$resolution_var,"_min_sample_",arguments$min_sample ,"_top_gene_cases.csv")), row.names = FALSE)
  return(cases_gene)
}

######Main######
dir.create(log_folder <- here("Data","model_case_sample_log"))
time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_")
logr::log_open(here(log_folder,paste0(time_case_prefix,"model_case_sample_logfile.log")))
logr::log_print(arguments)

logr::log_print(here("Results","CMH"))
xlsx <- here("Results","CMH",paste0("CMH_exact_summary_lclust_res_", arguments$resolution_var, "_min_sample_", arguments$min_sample, ".xlsx"))

cmhs <- xlsx %>%
    excel_sheets() %>%
    purrr::set_names() %>%
    map(read_excel, path = xlsx)

models <- names(cmhs) %>% print()

logr::log_print(models)

case_models <- lapply(models, get_cases)
names(case_models) <- models
write.xlsx(case_models, file = here("Results","CMH",paste0("cases_in_models_resolution_",arguments$resolution_var,"_min_sample_",arguments$min_sample ,".xlsx")))

logr::log_print("Finished")
logr::log_close()
