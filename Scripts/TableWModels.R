# TableWModels.R: Creates a table summarizing the models and parameters in analysis. Argument values cannot have spaces (for now, use _ instead of speace, will automaticaly  replace that with a space in the future). Values should be comma delimited.
'Usage: 
  TableWModels.R --model_name=<model_name> --inheritance=<inheritance> --effects=<effects>  [--case_group=<case_group>] [--footnote=<footnote>] [--gnomad_exome_pop_max=<gnomad_exome_pop_max>] [--gnomad_genome_pop_max=<gnomad_genome_pop_max>] [--mtr_max=<mtr_max>] [--revel_min=<revel_min>]
  
  Options:
  -h --help
    --model_name=<model_name> Example might be Ultra-Rare_Synonymous,Ultra-Rare_pLOF,Ultra-Rare_Missense,Ultra-Rare_Deleterious
    --inheritance=<inheritance> Dominant or recessive
    --effects=<effects> Examples might be Synonymous, pLOF, Missense
     --footnote=<footnote> example could be All pLOFs pass LOFTEE. All included variants have pext ratio > 0.9. See Methods for details. [default: NA]
     --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
     --gnomad_exome_pop_max=<gnomad_exome_pop_max> [default: blank]
     --gnomad_genome_pop_max=<gnomad_genome_pop_max> [default: blank]
     --mtr_max=<mtr_max> [default: blank]
     --revel_min=<revel_min> [default: blank]
' -> doc

library(docopt)
library(tidyverse)
library(data.table)
library(here)
library(logr)
library(stringr)
library(sjPlot)
"%!in%" <- Negate("%in%")
arguments <- docopt(doc, version = 'TableWModels.R 1.1')

# debug
# arguments <- list()
# arguments$model_name <- "Ultra-Rare_Synonymous,Ultra-Rare_pLOF,Ultra-Rare_Missense,Ultra-Rare_Deleterious"
# arguments$inheritance <- "Dominant,Dominant,Dominant,Dominant"
# arguments$effects <- "Synonymous,Missense,pLOF,Missense_and_pLOF"
# ************
# Setting up log and output folder
# ************
dir.create(here("Results"))
dir.create(output_folder <- here("Results","TableWModels"))
time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_",arguments$case_group, "_")
dir.create(here("Data"))
logr::log_open(here(paste0("Data/", time_case_prefix,"TableWModels.log")))
logr::log_print(arguments)
           
# ************
# Table with Models
# ************
write_df <- data.frame("Model Name" = strsplit(arguments$model_name, ",")[[1]],
                       "Inheritance" = strsplit(arguments$inheritance, ",")[[1]],
                       "Included Effects" = strsplit(arguments$effects, ",")[[1]],
                       check.names = FALSE)
if(arguments$gnomad_exome_pop_max != "blank") write_df$`gnomAD Exome Pop Max` <- strsplit(arguments$gnomad_exome_pop_max, ",")[[1]]
if(arguments$gnomad_genome_pop_max != "blank") write_df$`gnomAD Genome Pop Max` <- strsplit(arguments$gnomad_genome_pop_max, ",")[[1]]
if(arguments$mtr_max != "blank") write_df$`MTR` <- strsplit(arguments$mtr_max, ",")[[1]]
if(arguments$revel_min != "blank") write_df$`REVEL` <- strsplit(arguments$revel_min, ",")[[1]]
# names(write_df) <- str_replace(names(write_df), pattern = "\\.", replacement = " ")
library(sjPlot)

tab_df(write_df, file = paste0(output_folder, "/",time_case_prefix,"table_with_models.doc"), alternate.rows=TRUE, col.header = names(write_df),use.viewer = TRUE,
     show.footnote = TRUE, footnote = arguments$footnote)
# Close log
logr::log_close()


