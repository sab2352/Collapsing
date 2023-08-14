# regenieCovPheno.R performs three tasks. It extracts/reformats the principle components for the combined cluster calculated earlier. It creates a tab delimited interpretation of the vcf with just the ID and variant level annotations. TODO: can this be generated just from --list-var-anno with an ID list? I think so...
# Input: sample.txt, flashpca_eigenvectors, allQCBAFilterChr1variants.txt, allQCBAFilterChr2variants.txt, allQCBAFilterChr3variants.txt,allQCBAFilterChr4variants.txt, allQCBAFilterChr5variants.txt, allQCBAFilterChr6variants.txt, allQCBAFilterChr7variants.txt, allQCBAFilterChr8variants.txt, allQCBAFilterChr9variants.txt, allQCBAFilterChr10variants.txt, allQCBAFilterChr11variants.txt, allQCBAFilterChr12variants.txt, allQCBAFilterChr13variants.txt, allQCBAFilterChr14variants.txt, allQCBAFilterChr15variants.txt, allQCBAFilterChr16variants.txt, allQCBAFilterChr17variants.txt, allQCBAFilterChr18variants.txt, allQCBAFilterChr19variants.txt, allQCBAFilterChr20variants.txt, allQCBAFilterChr21variants.txt, allQCBAFilterChr22variants.txt, allQCBAFilterChrXvariants.txt
# Output: covariates.txt, phenotypes.txt, allQCBAFilterVariantMap.txt.gz, allQCBAFilterVariantMap_just_id.txt 

'Usage: 
  regenieCovPheno.R --pc_path=<pc_path> --coverage_path=<coverage_path> --regenie_folder=<regenie_folder> [--pcs=<pcs>] [--debug=<debug>]
  
  Options:
  -h --help
  --pc_path=<pc_path> principle component path
  --coverage_path=<coverage_path> coverage path
  --regenie_folder=<regenie_folder> folder with regenie project
  --pcs=<pcs> is the allele frequencies used for masking (includes singleton as default) [default: 6]
  --debug=<debug> [default: FALSE]
  #eventually will want the ability to do multiple phenotypes in the pheno file
  #best would be to use sample names to create phenotypes here - probably utilizing a new folder in /Regenie/
  
  
' -> doc

library(docopt)
library(tidyverse)
library(data.table)
library(broom)
library(here)
arguments <- docopt(doc, version = 'regenieCovPheno.R 1.1')

if(arguments$debug == "TRUE"){
  arguments <- list()
  arguments$pcs <- "6"
  arguments$pc_path <- here("Results","KinshipFlashPCA","PCAFlashLClust_res_0_4_cluster_11_min_sample_100_max_ratio_True_07_mixedCluster")
  arguments$coverage_path <- 
    here("Results","Coverage","CoverageFlashLClust_res_0_4_cluster_11_min_sample_100_max_ratio_True_07_mixedCluster")
  arguments$regenie_Folder <- here("Results","Regnie","0_4_cluster_11_min_sample_100_max_ratio_True")
}

######### Initiate Log ############
dir.create(log_folder <- here("Data","regenieCovPheno"))
time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_",arguments$case_group, "_")
logr::log_open(log_folder, paste0(time_case_prefix,"regenieCovPheno.log"))
logr::log_print(arguments)

####
# load in ped file from coverage folder
logr::log_print(sf <- list.files(here(arguments$coverage_path), "*_existing.sample.txt$", full.names = TRUE))
nrow(samples <- fread(sf, header = FALSE))
# load in PCs
logr::log_print(evf <- list.files(here(arguments$pc_path), "*_flashpca_eigenvectors$", full.names = TRUE))
nrow(pcs_prefiltered <- fread(evf))
# filter PCs for those in the ped file
nrow(pcs <- pcs_prefiltered %>% filter(FID %in% samples$V1))
# get pcs and sex as covariates
nrow(pcs_sex <- pcs %>% left_join(samples %>% select(V2, V5), by = c("IID" = "V2")) %>%
       rename(sex = V5))
max_pcs<-paste0("U",arguments$pcs)
nrow(cov_out <- pcs_sex %>% select(FID:all_of(max_pcs), sex))
# create phenotype file
pheno_all <- samples %>% rename(FID = V1, IID = V2) %>%
  mutate(All = V6 - 1) %>%
  select(FID, IID, All)
print(pheno_all %>% select(All) %>% table())
# save phenotype and covariate files
dir.create(new_dir <- here(arguments$regenie_folder,paste0("CovPheno_pc",arguments$pcs)))
write_delim(cov_out, here(new_dir, "covariates.txt"), delim = "\t")
write_delim(pheno_all, here(new_dir, "phenotypes.txt"), delim = "\t")
# save phenotype and covariate files
var.map.all <- c()

for (chr in c(1:22, "X")) {
  
  logr::log_print(paste0("chr ", chr))
  
  var.map <- fread(here(arguments$regenie_folder,"allQCBAFilter",paste0("allQCBAFilterChr", chr ,"variants.txt")))
  
  var.map.split <- var.map %>% separate(INFO, c("NS", "AF", "ANN"), sep = ";") %>%
    
    separate_rows(ANN, sep = ",") %>%
    
    separate(ANN, c("Effect", "Gene", "Transcript", "HGVS_c", "HGVS_p", "Polyphen_Humdiv", "Polyphen_Humvar"), sep = "\\|") %>%
    
    mutate(Effect = gsub("ANN=", "", Effect))
  
  
  
  var.map.all <- var.map.all %>% bind_rows(var.map.split)
  
}

write_delim(var.map.all, gzfile(here(arguments$regenie_folder,"allQCBAFilter","allQCBAFilterVariantMap.txt.gz")), delim = "\t")
write(x = var.map.all$ID ,  file = here(arguments$regenie_folder,"allQCBAFilter","allQCBAFilterVariantMap_just_id.txt"))

# create --aaf-file
# af=lapply(1:nrow(var.map.all), function(x) strsplit(var.map.all$AF[x], "="))

logr::log_print("Finished")
logr::log_close()

q()
