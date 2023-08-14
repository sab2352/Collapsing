# power_collapsing.R:  performs power or odds ratio calculations.

'Usage: 
  power_collapsing.R [--calc=<calc>] [--alpha=<alpha>] [--odds_ratio=<odds_ratio>] [--case_group=<case_group>] [--debug=<debug>]
  
  Options:
  -h --help
  --calc=<calc> what we are calculating, can be "power" or "es" (odds ratio) [default: power]
  --alpha=<alpha> type 1 error rate, or "Alpha" [default: 0.08]
  --odds_ratio=<odds_ratio> the detectable effect size (or odds ratio in the case of a binary outcome variable)[default: 3]
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
  --debug=<debug> debugger variable [default: FALSE]
  
' -> doc

library(docopt)
library(genpwr)
library(tidyverse)
library(data.table)
library(here)
library(logr)
arguments <- docopt(doc, version = 'power_collapsing.R 1.1')
#https://rdrr.io/cran/genpwr/f/vignettes/vignette.Rmd

if(arguments$debug == "TRUE"){
  arguments <- list()
  arguments$calc <- "power"
  arguments$alpha <- "0.0005"
  arguments$odds_ratio <- "11"
  arguments$case_group <- "temp"
  num_case <- 116
  num_total <- 10116
}

######### Initiate Log ############
log_folder <- here("Data","collapsing_powerLog")
if(!dir.exists(log_folder)) {
    dir.create(log_folder)
}
time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_",arguments$case_group, "_")
logr::log_open(here(log_folder,paste0( time_case_prefix,"collapsing_power_logfile.log")))
logr::log_print(arguments)
plot_path <- here("Results","Plots")
if(!dir.exists(plot_path)) {
    dir.create(plot_path)
}
logr::log_print(plot_path <- gsub("_/", "/", plot_path))

#get data
case_path <- list.files(here("Data"), "*case.seq.csv$", full.names = TRUE) 
full_data_path <- list.files(here("Data"), "*full_data.seq.csv$", full.names = TRUE)
logr::log_print(case_path)
logr::log_print(full_data_path)
num_case <- as.numeric(nrow(fread(case_path)))
logr::log_print(paste0("There are ", num_case, " cases"))
num_total <- as.numeric(nrow(fread(full_data_path)))
logr::log_print(paste0("There are ", num_total, " cases and controls"))

#produce graphs
if (arguments$calc == "power"){
  pwrare <- genpwr.calc(calc="power", model="logistic",ge.interaction=NULL, N = num_total, Case.Rate = (num_case / num_total),k=NULL, MAF = c(0.01,0.001,5E-4),OR = as.numeric(arguments$odds_ratio),Alpha=as.numeric(arguments$alpha), True.Model=c("Dominant", "Recessive", "Additive"), 
                        Test.Model=c("Dominant", "Recessive", "Additive", "2df"))
  pwcom <- genpwr.calc(calc="power", model="logistic",ge.interaction=NULL, N = num_total, Case.Rate = (num_case / num_total),k=NULL, MAF = c(0.1,0.05,0.01),OR = as.numeric(arguments$odds_ratio),Alpha=as.numeric(arguments$alpha), True.Model=c("Dominant", "Recessive", "Additive"), 
                       Test.Model=c("Dominant", "Recessive", "Additive", "2df"))
  power.plot(pwrare, y_log = T)
  ggsave(filename = here(plot_path, paste0(time_case_prefix, "power_rare.png")), device = "png", width = 10, height = 5)
  power.plot(pwcom, y_log = T)
  ggsave(filename = here(plot_path, paste0(time_case_prefix, "power_com.png")), device = "png", width = 10, height = 5)
}

if (arguments$calc == "es"){
  es <- genpwr.calc(calc="es", model="logistic",ge.interaction=NULL, N = num_total, Case.Rate = (num_case / num_total),k=NULL, MAF = c(0.1,0.05,0.01,0.001),Power = 0.8,Alpha=as.numeric(arguments$alpha), True.Model=c("Dominant", "Recessive", "Additive"), 
                    Test.Model=c("Dominant", "Recessive", "Additive", "2df"))
  or.plot(es,y_log = T)
  ggsave(filename = here(plot_path, paste0(time_case_prefix, "odds_ratio.png")), device = "png", width = 10, height = 5)
}

logr::log_print("Finished")
logr::log_close()
