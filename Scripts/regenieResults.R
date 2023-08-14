

'Usage: 
  regenieResults.R --regenie_folder=<regenie_folder> --pcs=<pcs> --vc_folder=<vc_folder> [--case_group=<case_group>] [--aafs=<aafs>] [--debug=<debug>]
  
  Options:
  -h --help
  --regenie_folder=<regenie_folder> folder with regenie project
  --pcs=<pcs>
  --vc_folder=<vc_folder>
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
  --aafs=<aafs> is the allele frequencies used for masking (includes singleton as default) [default: 0.1,0.05,0.01]
  --debug=<debug> [default: FALSE]
' -> doc
### to get qqman to work, need to run: /usr/local/igm/non-atav-tools/R-4.1.0_with_gcc_10-x86_64/bin/R then install.packages("qqman") to a temporary location in your folder.

library(docopt)
library(logr)
library(tidyr)
library(qqman)
library(here)
library(tidyverse)
library(data.table)
library(broom)
"%!in%" <- Negate("%in%")
# library(qqman,lib.loc="./R/x86_64-pc-linux-gnu-library/4.1/")
arguments <- docopt(doc, version = 'regenieResults.R 1.0')

if(arguments$debug == "TRUE"){
  arguments <- list()
  arguments$regenie_folder <- here("Results","Regenie","0_4_cluster_11_min_sample_100_max_ratio_True")
  arguments$vc_folder <- here(arguments$regenie_folder, "VC_pc6_mac_1")
  arguments$case_group <- "debug"
  arguments$pcs <- "6"
  arguments$aafs <- "0.1,0.05,0.01"
}



#########  ############
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
  initialize_logfile(arguments, "regenieResults")
}, error = function(e) {
  # Handle the error here
  print(paste("An error occurred while creating the logfile:", e$message))
})

#########  ############
aafs<-strsplit(arguments$aafs,split=",")[[1]]
aafs[length(aafs)+1]<-"singleton"
aafs[length(aafs)+1]<-"all"
masks<-fread(here(arguments$regenie_folder,"Anno","mask.mask"),select=1,header=F)
masks_full<-unlist(crossing(masks,aafs)%>%unite(.,col = "merged",sep = "."))

################################# Results ####################################

# packageurl <- "https://cran.r-project.org/src/contrib/Archive/GenABEL.data/GenABEL.data_1.0.0.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/GenABEL/GenABEL_1.8-0-6.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
# make sure you are in the right wd before loading here
library(GenABEL)


var.map.all <- fread(here(arguments$regenie_folder,"allQCBAFilter","allQCBAFilterVariantMap.txt.gz"))

####then pulls in all data for results####

for (dir2 in list.files(here(arguments$vc_folder))) {
  phenotypes <- fread(here(arguments$regenie_folder,paste0("CovPheno_pc", arguments$pcs),paste0("phenotypes",dir2,".txt")),nrows=0)
  phenotypes <- colnames(phenotypes)[3:length(colnames(phenotypes))]
  
  for (dir in list.files(here(arguments$vc_folder,dir2))) {
    for (pheno in phenotypes) {
      logr::log_print(sprintf("Reading from %s",dir))
      logr::log_print(sprintf("save_folder is %s",save_folder <- here(arguments$vc_folder,dir2,dir, paste0(pheno,"_",time_case_prefix,"output"))))
      dir.create(save_folder)
      logr::log_print(sprintf("Created directory %s",save_folder))
      reg_file <- grep("*.regenie",list.files(here(arguments$vc_folder,dir2,dir)),value=T)
      reg_file <- grep(pheno,reg_file,value = T)
      reg_file <-here(arguments$vc_folder,dir2,dir,reg_file)
      logr::log_print(sprintf("Reading from %s file",reg_file))
      res <- fread(reg_file,skip = 1)
      logr::log_print(sprintf("Finished reading %s",reg_file))
      if (nrow(res)==0){next}
      
      if (grepl("sum",dir)) {
        logr::log_print(sprintf("If statement is true"))
        res <- res %>% arrange(desc(LOG10P))
      
        res.filt <- res %>% filter(!is.na(LOG10P))
        
        res.filt <- res.filt %>% left_join(var.map.all %>% distinct(ID, Gene), by = c("ID" = "ID")) 
        
        res.filt %>% head()
      
        for (msk in masks_full){
          res.filt[which(grepl(msk,res.filt$ALLELE1)),]->res.filt.mask
          res.filt.mask$Pval<-10^(-(res.filt.mask$LOG10P))
          logr::log_print(sprintf("msk variable is %s", msk))
  
          reg_file_basename <- basename(reg_file)
          # write_csv(x = res.filt.mask, path = sprintf("%s%s",save_folder, gsub("\\.regenie", paste0("_",as.character(msk),".regenie"), reg_file_basename)))
          write_csv(x = res.filt.mask, file = sprintf("%s%s",save_folder,  gsub("\\.regenie", paste0("_",msk,".regenie"), reg_file_basename)))    
    
          logr::log_print(sprintf("Saving %s", filename_var  <- here(save_folder,paste0(dir2,"_" ,dir,"_",msk, "_manhattan.png"))))
          png(filename = filename_var, width = 800, height = 450)
          manhattan(res.filt.mask, chr="CHROM", bp="GENPOS" , snp="ID", p="Pval" )
          dev.off()
          
          lambda <- GenABEL::estlambda(res.filt.mask$Pval, method="regression")$estimate
          
          logr::log_print(sprintf("Saving %s", filename_var  <- here(save_folder,paste0(dir2,"_" ,dir,"_",msk, "_qq.png"))))
          qq(res.filt.mask$Pval, main = paste0("QQ Plot: Observed vs. expected p-values. Lambda = ", round(lambda, digits = 4)))
          dev.off()
        }
      } else {
        logr::log_print(sprintf("If statement is false"))
        res <- res %>% arrange(Pval)
      
        res.filt <- res %>% filter(!is.na(Pval))
        
        res.filt <- res.filt %>% left_join(var.map.all %>% distinct(ID, Gene), by = c("Name" = "ID")) 
        
        res.filt %>% head()
        
        ncases <- max(res.filt$Num_Cases, na.rm = TRUE)
        ncontrols <- max(res.filt$Num_Controls, na.rm = TRUE)
        
        # res.filt2 <- res.filt %>% filter(Num_Cases > ncases*0.9 & Num_Controls > ncontrols*0.9)
        # write_csv(res.filt2, gsub("\\.regenie", ".90pCov.regenie", reg_file))
        
        # res.filt3 <- res.filt2 %>% filter(Cases_Het + 2*Cases_Alt >= 3 & Controls_Het + 2*Controls_Alt >= 3)
        # write_csv(res.filt3, gsub("\\.regenie", ".90pCov.min_mac_case_control_3.regenie", reg_file))
        
        # 
        # png(filename = here("Regenie/Burden",dir2,dir,paste0(dir2,"_" ,dir,"_manhattan.png")), width = 800, height = 450)
        # manhattan(res.filt, chr="Chr", bp="Pos", snp="Name", p="Pval" )
        # dev.off()
        # 
        # lambda <- GenABEL::estlambda(res.filt$Pval, method="regression")$estimate
        # png(filename = here("Regenie/Burden",dir2,dir,paste0(dir2,"_" ,dir,"_qq.png")), width = 600, height = 600)
        # qq(res.filt$Pval, main = paste0("QQ Plot: Observed vs. expected p-values. Lambda = ", round(lambda, digits = 4)))
        # dev.off()
        # 
        # png(filename = here("Regenie/Burden",dir2,dir,paste0(dir2,"_" ,dir, "_90pCov_manhattan.png")), width = 800, height = 450)
        # manhattan(res.filt2, chr="Chr", bp="Pos", snp="Name", p="Pval" )
        # dev.off()
        # 
        # lambda <- GenABEL::estlambda(res.filt2$Pval, method="regression")$estimate
        # 
        # png(filename = here("Regenie/Burden",dir2,dir,paste0(dir2,"_" ,dir, "_90pCov_qq.png")), width = 600, height = 600)
        # qq(res.filt2$Pval, main = paste0("QQ Plot: Observed vs. expected p-values. Lambda = ", round(lambda, digits = 4)))
        # dev.off()
        # 
        # png(filename = here("Regenie/Burden",dir2,dir,paste0(dir2,"_" ,dir, "_90pCov.min_mac_case_control_3_manhattan.png")), width = 800, height = 450)
        # manhattan(res.filt3, chr="Chr", bp="Pos", snp="Name", p="Pval" )
        # dev.off()
        # 
        # lambda <- GenABEL::estlambda(res.filt3$Pval, method="regression")$estimate
        # 
        # png(filename = here("Regenie/Burden",dir2,dir,paste0(dir2,"_" ,dir, "_90pCov.min_mac_case_control_3_qq.png")), width = 600, height = 600)
        # qq(res.filt3$Pval, main = paste0("QQ Plot: Observed vs. expected p-values. Lambda = ", round(lambda, digits = 4)))
        # dev.off()
      reg_file_basename <- basename(reg_file)
        
      for (msk in masks_full){
        ##can change res.filt2 to either res.filt2 or res.filt3 if you want to get rid of the coverage harmonization of 90% or add a min Case/control number
        # the code filters the rows of res.filt where the Alt column contains the string msk and assigns the filtered subset to a new object called res.filt.mask.
        logr::log_print(sprintf("msk variable is %s", msk))
        res.filt[which(grepl(msk,res.filt$Alt)),]->res.filt.mask.allModels
        for (model in unique(res.filt.mask.allModels$Model)){
          logr::log_print(sprintf("model variable is %s", model))
          res.filt.mask <- filter(res.filt.mask.allModels, Model == model)
          
          base_save <- here(save_folder, paste0(reg_file_basename, "_",msk,"_",model))
          write_csv(x = res.filt.mask, file = paste0(base_save, ".csv"))
          # write_csv(x = res.filt.mask, file = sprintf("%s%s",save_folder, gsub("\\.regenie", paste0("_",as.character(msk),".regenie"), reg_file_basename)))
          logr::log_print(sprintf("Saving %s",filename_var <- paste0(base_save, "_manhattan.png")))
          png(filename = filename_var, width = 800, height = 450)
          manhattan(res.filt.mask, chr="Chr", bp="Pos", snp="Name", p="Pval" )
          dev.off()
          
          lambda <- GenABEL::estlambda(res.filt.mask$Pval, method="regression")$estimate
          
          logr::log_print(sprintf("Saving %s", filename_var <- paste0(base_save, "_qq.png")))
          png(filename = filename_var, width = 600, height = 600)
          qq(res.filt.mask$Pval, main = paste0("QQ Plot: Observed vs. expected p-values. Lambda = ", round(lambda, digits = 4)))
          dev.off()
        }
       }
      }
    }
  }
}

