# cmh_permute_lclust_short_script.R: runs permutations. Created by Gundula. Adapted by Josh and others.
library(optparse)
library(logr)
library(here)
library(docopt)

GetOptList <- function() {
  # Get option list.
  #
  # Returns:
  #   Option list to parse options.
  option.list <- list(
    make_option(c("-p", "--project"), action="store", 
                default="", dest="project", 
                help="path to Results directory of project, [default %default]"),

    make_option(c("-d", "--dir"), action="store", 
                default="", dest="dir", 
                help="name of group, [default %default]"),

    make_option(c("-r", "--res"), action="store", 
                default="0_3", dest="res", 
                help="resolution, [default %default]"),

    make_option(c("-n", "--minsample"), action="store", 
                default=5, dest="minsample", 
                help="minimum number of cases/controls per cluster, [default %default]"),

    make_option(c("-m", "--model"), action="store", 
                default="dominantSynonymous", dest="model", 
                help="model to permute, [default %default]"),

    make_option(c("-s", "--permstart"), action="store", 
                default=1, dest="permstart", 
                help="first permutation, [default %default]"),

    make_option(c("-e", "--permend"), action="store", 
                default=1, dest="permend", 
                help="last permutation, [default %default]"),

    make_option(c("-c", "--cores"), action="store", 
                default=5, dest="cores", 
                help="number of cores, [default %default]")

  )
}

option.list <- GetOptList()
opt <- parse_args(OptionParser(option_list=option.list))

# debug
# opt <- list()
# opt$cores <- 1
# opt$permend <- 50
# opt$permstart <- 1
# opt$model <- "URFUNC_min80_igm"
# opt$minsample <- 20
# opt$res <- "0_4"
# opt$dir <- ""

source(here("Scripts","cmh_permute_lclust_short_functions.R"))
"%!in%" <- Negate("%in%")

dir.create(here("Data"))
dir.create(log_folder <- here("Data","cmh_permute_lclust_short_script_log"))
time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_")
logr::log_open(here(log_folder, paste0(time_case_prefix,"cmh_permute_lclust_short_scripts_logfile.log")))
logr::log_print(opt)
logr::log_print(cores <- opt$cores)

logr::log_print(sprintf("Analyzing the following models %s", model_list <- unlist(strsplit(opt$model,","))))
logr::log_print(dir <- opt$dir)



if (dir == "") {
    dir = ""
} else {
    dir <- paste0(dir, "_")
}


pc_path <- here("Results","KinshipFlashPCA", dir)
plot_path <- here(paste0("Results/Plots_", gsub("_", "",  dir), "/"))
plot_path <- gsub("_/", "/", plot_path)
results_path <- here(paste0("Results/CMH_", gsub("_", "",  dir), "/"))
results_path <- gsub("_/", "/", results_path)
permut_path <- here(paste0("Results/Permutations_", gsub("_", "",  dir), "/"))
permut_path <- gsub("_/", "/", permut_path)


if(!dir.exists(permut_path)) {
  dir.create(permut_path)
} else {
  logr::log_print((paste0("Directory ", permut_path, "already exists")))
}

if(!dir.exists(here("Data","cmh_permutation_lclust_Log"))) {
  dir.create(here("Data","cmh_permutation_lclust_Log"))
} else {
  logr::log_print((paste0("Directory Log already exists")))
}

resolution <- opt$res
min_sample <- opt$minsample

# adapt to your filename
logr::log_print(cl_size_path <- list.files(plot_path, pattern = paste0("*lclustering_res_", resolution, "_cluster_sizes.txt$"), full.names = TRUE))
cl_sizes <- fread(cl_size_path)
logr::log_print("Analyzing the following clusters")
logr::log_print(nclust <- (cl_sizes %>% filter(case >= min_sample & control >= min_sample))$cluster)

# you can start with 1:100 and later run 101:1000 if you want since the number is used as seed
model <- model_list[1]
for (model in model_list){
  logr::log_print(sprintf("Now analyzing %s", model))
  
  for(i in 1:2){
    nperm <- opt$permstart:opt$permend
    files <- list.files(permut_path, pattern = glob2rx(paste0(dir, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_", "*.RDS")))
    success <- as.numeric(gsub(paste0(dir, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_|.RDS"), "", files))
    logr::log_print(head((nperm <- nperm[!nperm%in%success])))
    if(length(nperm) > 0){
      logr::log_print(sprintf("Now running getCMHresP on model %s", model))
      # run for each model separately
      getCMHresP(model, dir, nclust, nperm, cores)
    }
  }
}
# close log file
logr::log_close()

# getCMHresP(model, dir, pc_path, nclust, nperm, cores)

# to check if there are permutations missing
# files <- list.files(paste0(permut_path), pattern = glob2rx(paste0(dir, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_", "*.RDS")))
# success <- as.numeric(gsub(paste0(dir, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_|.RDS"), "", files))
# nperm <- nperm[!nperm%in%success]
