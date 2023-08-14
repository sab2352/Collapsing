# cmh_permute_lclust_short.R: runs permutations. Created by Gundula. Adapted for atav by Josh and others.
# THIS IS UNDER CONSTRUCTION. STARTING WITH MODIFICATIONS OF SCRIPT FIRST TO TAKE ADVANTAGE OF QSUB

'Usage: 
  cmh_permute_lclust_short.R --resolution_var=<resolution_var> --min_sample=<min_sample> [--models_to_exclude=<models_to_exclude>]
  
  Options:
  -h --help
  --min_sample=<min_sample> minimum number of cases/controls of included clusters
  --resolution_var=<resolution_var> resolution for clustering.
  --models_to_exclude=<models_to_exclude> exclude these models. Typically the master genotype model. [default: dominantNoneMAF]
  
' -> doc

library(here)
library(docopt)
library(data.table)
source(here("Scripts/cmh_permute_lclust_short_functions.R"))
"%!in%" <- Negate("%in%")

arguments <- docopt(doc, version = 'cmh_permute_lclust_short.R 1.1')


dir <- ""

pc_path <- paste0("KinshipFlashPCA", dir)

if (dir == "") {
    dir = ""
} else {
    dir <- paste0(dir, "_")
}

plot_path <- here(paste0("Results/Plots_", gsub("_", "",  dir), "/"))
plot_path <- gsub("_/", "/", plot_path)
results_path <- here(paste0("Results/CMH_", gsub("_", "",  dir), "/"))
results_path <- gsub("_/", "/", results_path)
permut_path <- here(paste0("Results/Permutations_", gsub("_", "",  dir), "/"))
permut_path <- gsub("_/", "/", permut_path)


if(!dir.exists(permut_path)) {
  dir.create(permut_path)
} else {
  print((paste0("Directory ", permut_path, "already exists")))
}

if(!dir.exists("Log")) {
  dir.create("Log")
} else {
  print((paste0("Directory Log already exists")))
}

resolution <- arguments$resolution_var
min_sample <- as.numeric(arguments$min_sample)

# adapt to your filename
print(cl_size_path <- list.files(plot_path, pattern = paste0("*lclustering_res_", resolution, "_cluster_sizes.txt$"), full.names = TRUE))
cl_sizes <- fread(cl_size_path)
nclust <- (cl_sizes %>% filter(case >= min_sample & control >= min_sample))$cluster %>% print()

model <- "dominantPTV"
model <- "dominantSynonymous"
model <- "dominantUltraRareP"

model <- "dominantRareNB"
model <- "dominantRareEnsemble"
model <- "dominantUltraRareEnsemble"

# you can start with 1:100 and later run 101:1000 if you want since the number is used as seed
nperm <- 1:250
cores <- 5

# run for each model separately
getCMHresP(model, dir, pc_path, nclust, nperm, cores)

# to check if there are permutations missing
files <- list.files(paste0(permut_path), pattern = glob2rx(paste0(dir, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_", "*.RDS")))
success <- as.numeric(gsub(paste0(dir, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_|.RDS"), "", files))
nperm <- nperm[!nperm%in%success]
