library(optparse)

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


setwd(opt$project)
source("../Scripts/cmh_permute_lclust_short_genderstrat_functions.R")

dir <- opt$dir

pc_path <- paste0("KinshipFlashPCA", dir)

if (dir == "") {
    dir = ""
} else {
    dir <- paste0(dir, "_")
}


plot_path <- paste0("Plots_", gsub("_", "",  dir), "/")
plot_path <- gsub("_/", "/", plot_path)
results_path <- paste0("CMH_", gsub("_", "",  dir), "/")
results_path <- gsub("_/", "/", results_path)
permut_path <- paste0("Permutations_", gsub("_", "",  dir), "/")
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

resolution <- opt$res
min_sample <- opt$minsample
genders <- c("male", "female")

# adapt to your filename
cl_sizes <- fread(paste0("../", plot_path, dir, "lclustering_res_", resolution, "_cluster_sizes_gender.txt"))

nclusts <- c()
for (gender in genders) {
  nclusts[[gender]] <- (cl_sizes %>% 
                        filter(Gender == gender & case >= min_sample & control >= min_sample))$cluster
}

model <- opt$model

# model <- "dominantPTV"
# model <- "dominantSynonymous"
# model <- "dominantUltraRareP"

# model <- "dominantRareNB"
# model <- "dominantRareEnsemble"
# model <- "dominantUltraRareEnsemble"

# you can start with 1:100 and later run 101:1000 if you want since the number is used as seed
nperm <- opt$permstart:opt$permend
cores <- opt$cores

# run for each model separately
getCMHresP(model, dir, genders, pc_path, nclusts, nperm, cores)

# to check if there are permutations missing
# files <- list.files(paste0(permut_path), pattern = glob2rx(paste0(dir, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_", "*.RDS")))
# success <- as.numeric(gsub(paste0(dir, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_|.RDS"), "", files))
# nperm <- nperm[!nperm%in%success]
