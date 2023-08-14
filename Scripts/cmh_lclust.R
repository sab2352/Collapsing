# cmh_lclust.R: runs CMH on clusters. Created by Gundula. Adapted for atav by Josh and others.
# make sure you are in the right wd before loading here

'Usage: 
  cmh_lclust.R --resolution_var=<resolution_var> --min_sample=<min_sample> [--models_to_exclude=<models_to_exclude>] [--debug=<debug>]
  
  Options:
  -h --help
  --min_sample=<min_sample> minimum number of cases/controls of included clusters
  --resolution_var=<resolution_var> resolution for clustering.
  --models_to_exclude=<models_to_exclude> exclude these models. Typically the master genotype model. [default: dominantNoneMAF]
  --debug=<debug> [default: FALSE]
' -> doc

library(here)
library(docopt)
library(data.table)
library(logr)
source(here("Scripts","cmh_lclust_functions.R"))
source(here("Scripts","getGeneGenos_functions.R"))
"%!in%" <- Negate("%in%")

arguments <- docopt(doc, version = 'cmh_lclust.R 1.1')
if(arguments$debug == "TRUE"){
  arguments <- list()
  arguments$min_sample <- "15"
  arguments$resolution_var <- "0_4"
  arguments$models_to_exclude <- "dominantNoneMAF,recessiveAutosomalwoMAPIDHP,recessiveAutosomalLEwoMAPIDHP"
}

#########  ############
dir.create(log_folder <- here("Data","cmh_lclust_Log"))
time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_")
logr::log_open(paste0(log_folder,"/", time_case_prefix,"cmh_lclust_logfile.log"))
logr::log_print(arguments)
logr::log_print(sprintf("models to exclude %s",models_to_exclude <- unlist(strsplit(arguments$models_to_exclude,","))))

#########  ############
resolution <- arguments$resolution_var
min_sample <- as.numeric(arguments$min_sample)

dir <- ""

if (dir == "") {
    dir = ""
} else {
    dir <- paste0(dir, "_")
}

plot_path <- here(paste0("Results/Plots_", gsub("_", "",  dir), "/"))
plot_path <- gsub("_/", "/", plot_path)
results_path <- here(paste0("Results/CMH_", gsub("_", "",  dir), "/"))
results_path <- gsub("_/", "/", results_path)


if(!dir.exists(results_path)) {
    dir.create(results_path)
} else {
  logr::log_print((paste0("Directory ", results_path, "already exists")))
}

# load in cluster sizes
logr::log_print(cl_size_path <- list.files(plot_path, pattern = paste0("*lclustering_res_", resolution, "_cluster_sizes.txt$"), full.names = TRUE))
cl_sizes <- fread(cl_size_path)
nclust <- (cl_sizes %>% filter(case >= min_sample & control >= min_sample))$cluster %>% print()
# nclust <- c(0, 1, 3) # if you want to select your own clusters to combine

# model <- "dominantUltraRareEnsemble"
# cmh <- getCMHres(model, dir, nclust) # in case you want to run a single model

logr::log_print("Folders of included clusters")
logr::log_print(nclust_folders <- lapply(nclust, function(cl) here("Results", paste0("Collapsing", gsub("_", "", dir)), paste0("LClust_res_", resolution, "_cluster_", cl, "_", dir, "FlashColl_07"))))
logr::log_print("Models in folders")
logr::log_print(models <- lapply(nclust_folders, function(cl_folder) list.dirs(cl_folder, recursive = FALSE, full.names = FALSE)))
logr::log_print("Common models")
logr::log_print(models.r <- Reduce(intersect, models))
logr::log_print("After excluding models")
logr::log_print(models.r <- models.r[models.r %!in% models_to_exclude])

# logr::log_print(sprintf("Analyzing the following models %s", models.r <- models.r[-grep("None", models.r)]))
# logr::log_print(sprintf("Analyzing the following models %s", models.r <- models.r[-grep("None|recessive", models.r)]))

logr::log_print(sprintf("Models to create %s", models_to_create <- paste0(models.r,"_res_", resolution, "_cluster_", paste(nclust, collapse = "_"), "_exact_CMH.RDS")))
logr::log_print(sprintf("Models already completed %s",models_already_completed <- list.files(here("Results","CMH"), pattern = paste0("_res_", resolution, "_cluster_", paste(nclust, collapse = "_"), "_exact_CMH.RDS$"))))
logr::log_print("Models which will run to create the .RDS files")
logr::log_print(models_to_run_paths <- models_to_create[models_to_create %!in% c(models_already_completed)])
logr::log_print("Extracting model names")
logr::log_print(models_to_run_list <- strsplit(models_to_run_paths, paste0("_res_", resolution)))
if(length(models_to_run_list) > 0){ 
  logr::log_print(sprintf("Models to run %s",models_to_run <- sapply(1:length(models_to_run_list), function(x) models_to_run_list[[x]][1])))
  cmh.all <- lapply(models_to_run[models_to_run %!in% models_to_exclude], function(x) getCMHres(x, dir, nclust))
}


loadCMH <- function(model, dir) {
  print(model)
  cmh <- readRDS(paste0(results_path, "/", dir, model, "_res_", resolution, "_cluster_", paste(nclust, collapse = "_"), "_exact_CMH.RDS"))
  return(cmh)
}

cmh.all <- lapply(models.r[models.r %!in% models_to_exclude], function(x) loadCMH(x, dir))


names(cmh.all) <- gsub("dominant", "" , models.r)
logr::log_print("Writing excel")
write.xlsx(cmh.all, here(results_path, paste0("CMH_exact_summary_lclust_res_", resolution, "_min_sample_", min_sample, ".xlsx")))


# close log file
logr::log_close()
stop("Finished")
# Have not debugged beyond here

# all models in one table
cmh <- unnest(tibble(cmh.all), .id = "Model")

# top 100 genes over all models
cmh %>% filter(!str_detect(Model, "Synonymous")) %>% 
    arrange((p.value)) %>% distinct(gene, .keep_all = TRUE) %>% 
    select(Model, gene:`%QV+ Ctrl`) %>% top_n(-100, wt = p.value) %>% as.data.frame()



# for extracting a single gene from the genotypes files of all clusters (original function also works with gene list)
gene <- "NPHS2"

model <- "dominantUltraRareEnsemble"

geneGeno <- getGeneListGenos(model, dir, nclust, gene, recessive = FALSE)[[1]] %>% 
    mutate_at(vars(contains("HGMD|ClinVar")), as.character)
write_csv(geneGeno, paste0(results_path, "/", dir, model, "_genotypes_", gene, "_res_", resolution, "_min_sample_", min_sample, ".csv"))

model <- "recessiveAutosomalLEwoMAPIDHP"
geneGeno <- getGeneListGenos(model, dir, nclust, gene, recessive = TRUE)[[1]] %>% 
    mutate_at(vars(contains("HGMD|ClinVar")), as.character)

write_csv(geneGeno, paste0(results_path, "/", dir, model, "_genotypes_", gene, "_res_", resolution, "_min_sample_", min_sample, ".csv"))
