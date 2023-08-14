# make sure you are in the right wd before loading here
library(here)

source(here("Scripts/cmh_lclust_gender_functions.R"))

resolution <- "0_3"
min_sample <- 5

dir <- "All"

if (dir == "") {
    dir = ""
} else {
    dir <- paste0(dir, "_")
}

plot_path <- here(paste0("Plots_", gsub("_", "",  dir), "/"))
plot_path <- gsub("_/", "/", plot_path)
results_path <- here(paste0("Results/CMH_", gsub("_", "",  dir), "/"))
results_path <- gsub("_/", "/", results_path)

if(!dir.exists(results_path)) {
    dir.create(results_path)
} else {
    print((paste0("Directory ", results_path, " already exists")))
}


# adapt path
  cl_sizes <- fread(paste0(plot_path, "/", dir, "lclustering_res_", resolution, "_cluster_sizes_gender.txt"))

gender <- "male"
nclust <- (cl_sizes %>% filter(Gender == gender & case >= min_sample & control >= min_sample))$cluster %>% print()
# nclust <- c(0, 1, 3) # if you want to select your own clusters to combine

# model <- "dominantUltraRareEnsemble"
# cmh <- getCMHres(model, dir, nclust, gender) # in case you want to run a single model

# adapt path
models <- lapply(nclust, function(cl) list.dirs(here(paste0("Results/Collapsing", gsub("_", "", dir), "/LClust_res_", resolution, "_cluster_", cl, "_", dir, "FlashColl_07_", gender)), recursive = FALSE, full.names = FALSE))
models.r <- Reduce(intersect, models)
# models.r <- models.r[-grep("None|recessive", models.r)]
models.r <- models.r[grep("recessive", models.r)]

cmh.all <- lapply(models.r, function(x) getCMHres(x, dir, nclust, gender))
names(cmh.all) <- gsub("dominant", "" , models.r)
write.xlsx(cmh.all, paste0(results_path, "/", dir, "CMH_exact_summary_lclust_res_", resolution, "_min_sample_", min_sample, "_", gender, "_recessive.xlsx"))

gender <- "female"
nclust <- (cl_sizes %>% filter(Gender == gender & case >= min_sample & control >= min_sample))$cluster %>% print()

cmh.all <- lapply(models.r, function(x) getCMHres(x, dir, nclust, gender))
names(cmh.all) <- gsub("dominant", "" , models.r)
write.xlsx(cmh.all, paste0(results_path, "/", dir, "CMH_exact_summary_lclust_res_", resolution, "_min_sample_", min_sample, "_", gender, "_recessive.xlsx"))

genders <- c("male", "female")
nclusts <- list()
for (gender in genders) {
  nclusts[[gender]] <- (cl_sizes %>% 
                        filter(Gender == gender & case >= min_sample  & control >= min_sample))$cluster
}

cmh.all <- lapply(models.r, function(model) getCMHresGenderStrat(model, dir, nclusts, genders))
names(cmh.all) <- gsub("dominant", "", models.r)

write.xlsx(cmh.all, paste0(results_path, "/", dir, "CMH_exact_summary_lclust_res_", resolution, "_min_sample_", min_sample, "_genderstrat_recessive.xlsx"))

