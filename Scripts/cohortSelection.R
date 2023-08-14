# cohortSelection.R: creates a cohort, Created by Gundula. Adapted for atav by Josh and others.

'Usage: 

  cohortSelection.R --case_list_path=<case_list_path> [--ctrl_list_path=<ctrl_list_path>] [--sample_path_w_filename=<sample_path_w_filename>] [--case_group=<case_group>] [--exclude_list=<exclude_list>] [--BroadPhenotype_include=<BroadPhenotype_include>] [--acceptable_capture_kits=<acceptable_capture_kits>] [--DetailedPhenotype_excl_regex=<DetailedPhenotype_excl_regex>] [--CCDSBasesCov10X_var=<CCDSBasesCov10X_var>] [--DBSNPOverlapSNVs_var=<DBSNPOverlapSNVs_var>] [--DBSNPOverlapIndels_var=<DBSNPOverlapIndels_var>] [--ContaminationPercent_var=<ContaminationPercent_var>] [--Ancestry_var=<Ancestry_var>] [--Ancestry_filter=<Ancestry_filter>]

  
  Options:
  -h --help
  --case_list_path=<case_list_path> name of text file without header with list of cases to be included. Must be located in Input folder at top level of repo
  --ctrl_list_path=<ctrl_list_path> path to text file without header with list of ctrls to be included [default: NA]  
  --sample_path_w_filename=<sample_path_w_filename> path to csv file downloaded from https://sequence.igm.cumc.columbia.edu/ [default: /nfs/goldstein/software/atav_home/data/sample/igm_lims_sample_master_list.tsv]
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
  --exclude_list=<exclude_list> specific samples to exclude [default: NA]  
  --BroadPhenotype_include=<BroadPhenotype_include> broad phenotypes to include in controls. Must uses _ instead of space. [default: healthy_family_member,control,control_mild_neuropsychiatric_disease]
  --acceptable_capture_kits=<acceptable_capture_kits> acceptable capture kits [default: NexteraRapidCapture,65MB,Roche,RocheV2,IDTERPv1,IDTERPv1Plus,IDTxGEN,IDTERPv1mtDNA,AgilentCRE,AgilentV4,Agilentv5,AgilentV5,AgilentV5UTR,AgilentV6,MedExome,N/A,Genome_v1,IDTERPv2]
  --DetailedPhenotype_excl_regex=<DetailedPhenotype_excl_regex> regexp string to exclude things in detailed phenotype [default: NA] 
  --CCDSBasesCov10X_var=<CCDSBasesCov10X_var> [default: 90] samples must have at least this % of the CCDS region covered by 10x
  --DBSNPOverlapSNVs_var=<DBSNPOverlapSNVs_var> [default: 0.85] must have at least this % overlap of SNVs with DBSNP
  --DBSNPOverlapIndels_var=<DBSNPOverlapIndels_var> [default: 0.80] must have at least this % overlap of indels with DBSNP
  --ContaminationPercent_var=<ContaminationPercent_var> [default: 2] upper limit to the contamination percentage for samples to include. Ignored if coverage >200x at 15% is used for these
  --debug=<debug> [default: FALSE]
  --Ancestry_filter=<Ancestry_filter> [default: NA] Predicted ancestry you wish to filter on
  --Ancestry_var=<Ancestry_var> [default: 0.75] Min ancestry percentage 

' -> doc
library(docopt)
library(tidyverse)
library(data.table)
library(readxl)
library(here)
library(logr)
"%!in%" <- Negate("%in%")
arguments <- docopt(doc, version = 'cohortSelection 1.1')
# debugging
if(arguments$debug){
  arguments <- list()
  arguments$exclude_list <- "IGM-GHARMaGICP22uuA038fSF068,IGM-GHARMaGICP22uuA151fSF231,IGM-GHARMaGIC22uuSF036-II-1-P2fSF036"
# arguments$case_list_path <- "cohort_in_dragen.txt"
# arguments$exclude_list <- "NA"
# arguments$ctrl_list_path <- "NA"
# arguments$acceptable_capture_kits <- "65MB,Roche,RocheV2,IDTERPv1,IDTERPv1Plus,IDTxGEN,IDTERPv1mtDNA,AgilentCRE,AgilentV4,Agilentv5,AgilentV5,AgilentV5UTR,AgilentV6,MedExome,N/A,Genome_v1"
# arguments$BroadPhenotype_include <- "healthy_family_member,control,kidney_and_urological_disease,amyotrophic_lateral_sclerosis,obsessive_compulsive_disorder,dementia,epilepsy,liver_disease,other_neurological disease,other_neuropsychiatric_disease,other_neurodevelopmental_disease,intellectual_disability,control_mild_neuropsychiatric_disease"
# arguments$DetailedPhenotype_excl_regex <- "epilep|seizure"
# arguments$sample_path_w_filename<-"/Users/jm4279/OneDrive - cumc.columbia.edu/sequence_files/20221001_indragendb.csv"
# arguments$ContaminationPercent_var <- "2"
# arguments$CCDSBasesCov10X_var <- "90"
# arguments$DBSNPOverlapSNVs_var <- "0.85"
# arguments$DBSNPOverlapIndels_var <- "0.80"
# arguments$Ancestry_filter <- "NA"
# arguments$Ancestry_var <- "0.75"
}


# Create log file and cast variables
tryCatch({
  dir.create(here("Data"))
  dir.create(log_folder <- here("Data","cohortSelectionLog"))
  time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_",arguments$case_group, "_")
  logr::log_open(here(log_folder,paste0( time_case_prefix,"cohortSelection_logfile.log")))
  logr::log_print(arguments)
  logr::log_print("Case list...")
}, error = function(e) {
  message("Caught an error while creating log file: ", e$message)
})

tryCatch({
  if (ncol(case_list <- fread(here("Input",arguments$case_list_path), header = FALSE)) == 8) {
    inputfam <- case_list
    logr::log_print(case_list <- subset(inputfam, V6 == 2))
    ctrl_list <- subset(inputfam, V6 == 1)
    fam_file <- TRUE
  } else {
    logr::log_print(case_list)
    fam_file <- FALSE
  }
}, error = function(e) {
  message("Caught an error while loading in case list: ", e$message)
})

tryCatch({
  logr::log_print(acceptable_capture_kits <- unlist(strsplit(arguments$acceptable_capture_kits,",")))
  logr::log_print(BroadPhenotype_include <- unlist(strsplit(gsub("_"," ", arguments$BroadPhenotype_include),",")))
  logr::log_print(CCDSBasesCov10X_var <- as.numeric(arguments$CCDSBasesCov10X_var))
  logr::log_print(DBSNPOverlapSNVs_var <- as.numeric(arguments$DBSNPOverlapSNVs_var))
  logr::log_print(DBSNPOverlapIndels_var <- as.numeric(arguments$DBSNPOverlapIndels_var))
  logr::log_print(ContaminationPercent_var <- as.numeric(arguments$ContaminationPercent_var))
  logr::log_print(exclude_list_additions <- unlist(strsplit(arguments$exclude_list,",")))
  
  exclude_list <- c("CGNDHDA01633", "CGNDHDA01594","Diagseq3531f1493","washei49207", "washei49209", "washei49245", "washei49300",
                    "washei50375", "washei50449", "washei50459", "washei50470", "washei50476", "washei50477", "washei50494",
                    "GHARRHDY18uuUROGE313xx2","GHARRHDY18uuUROGE434xx3", "GHARRHDY18uuP1027xx3", "GHARRHDY18uuUROME108xx2",
                    "GHARRHDY18uuUROGE211xx1", "GHARRHDY18uuUROGE225xx2", "GHARRHDY18uuUROGE246xx2", "GHARRHDY18uuUROGE138xx1",
                    "fetal0119M", "fetal0119D", "fetal0120M", "fetal0120D", "fetal0116M", "fetal0118M", "fetal0118D", "fetal0125D",
                    "fetal0127D", "fetal0122M", "Diagseq3330f1341", "sqcudn725547",exclude_list_additions) #defaults set by gundula and added by josh. 
  logr::log_print("Diagseq3531f1493, has a problem with coverage bins in dragen. see ticket https://redmine.igm.cumc.columbia.edu/issues/7142")
  logr::log_print(paste0("Exclude list ", exclude_list))
}, error = function(e) {
  message("Caught an error casting variables: ", e$message)
})

tryCatch({
  ######### change this ############
  # sample_path <- "pathtosamples"
  sample_path <- dirname(arguments$sample_path_w_filename)
  logr::log_print(nrow(data <- fread(arguments$sample_path_w_filename))) # csv file downloaded from https://sequence.igm.cumc.columbia.edu/
  
  logr::log_print(paste0("Number of samples in unfiltered sequence data: ", nrow(data %>% distinct(sample_internal_name))))
  
  cases_not_in_dragen <- case_list %>% filter(tolower(V1) %!in% tolower(data$sample_internal_name))
  if(nrow(cases_not_in_dragen)){
    logr::log_print(paste0("There are cases absent from dragen: ", cases_not_in_dragen$V1))
    write.table(file = here("Data", paste0(time_case_prefix,"cases_absent_from_dragen.txt")), 
                x = cases_not_in_dragen$V1, row.names = FALSE, col.names = FALSE, quote = FALSE) 
  }
  
  data <- data %>% 
      mutate(capture_kit =  case_when(capture_kit == "" | is.na(capture_kit) | is.null(capture_kit) ~ exomeKit, TRUE ~ capture_kit)) #capture kit is generated by pipeline, exomeKit is submitted by user. Capture_kit is used to initialize "atavdb capture kit" which is blind to user but part of atav processing
  
  # following step can be skipped since we remove them in kinship step, but to get a better idea of final control numbers we already remove gnomAD samples here
  gnomAD <- fread(here("DefaultData","igm_gnomad_sample_062620.txt"), header = FALSE) # copied from /nfs/goldstein/software/atav_home/data/sample/igm_gnomad_sample_062620.txt
  
  filtData <- data %>% 
      filter(sample_status == "In DragenDB") %>% 
      filter(sample_type %in% c("Exome", "Genome")) %>%
      filter(tolower(sample_internal_name) %!in% tolower(exclude_list)) %>% # there is something wrong with these 2 samples
      filter(ContaminationPercentage <= ContaminationPercent_var | (ContaminationPercentage <= 15 & MeanCoverage > 200)) %>% # can be relaxed. Per Gundula "The contamination calculation does not work for high coverage with bwa-mem type aligners according to Dan. If you check the lane information in RG metrics. These samples are not contaminated."
      filter(CCDSBasesCov10X >= CCDSBasesCov10X_var) %>% # # can be relaxed 
      filter(tolower(capture_kit) %in% tolower(acceptable_capture_kits)) %>% # add additional if you have cases that have other values and want to include them
      filter(DBSNPOverlapSNVs >= DBSNPOverlapSNVs_var & DBSNPOverlapIndels >= DBSNPOverlapIndels_var) %>% # can be relaxed
      filter(!sample_internal_name %in% gnomAD$V1) %>%
      identity()    
  
  if (arguments$Ancestry_filter != "NA") {
    logr::log_print(paste0("Limiting Ancestry to ", arguments$Ancestry_filter))
    if (arguments$Ancestry_filter == "European") {
      filtData <- filtData %>% filter(Caucasian_prob > arguments$Ancestry_var & genotyping_rate != -1)
    } else if (arguments$Ancestry_filter == "MiddleEastern") {
      filtData <- filtData %>% filter(MiddleEastern_prob > arguments$Ancestry_var & genotyping_rate != -1)
    } else if (arguments$Ancestry_filter == "Latino") {
      filtData <- filtData %>% filter(Hispanic_prob > arguments$Ancestry_var & genotyping_rate != -1)
    } else if (arguments$Ancestry_filter == "EastAsian") {
      filtData <- filtData %>% filter(EastAsian_prob > arguments$Ancestry_var & genotyping_rate != -1)
    } else if (arguments$Ancestry_filter == "SouthAsian") {
      filtData <- filtData %>% filter(SouthAsian_prob > arguments$Ancestry_var & genotyping_rate != -1)
    } else if (arguments$Ancestry_filter == "African") {
      filtData <- filtData %>% filter(African_prob > arguments$Ancestry_var & genotyping_rate != -1)
    } else {
      logr::log_print("Ancestry_filter does not match any ancestry")
    }
  }

  logr::log_print(paste0("Number of samples in filtered dragen db: ", nrow(filtData %>% distinct(sample_internal_name))))
  
  # adapt based on your phenotype
  controls_1 <- filtData %>% 
    filter(tolower(sample_internal_name) %!in% tolower(case_list$V1), toupper(seqGender) %in% c("M", "F")) %>%
    filter(tolower(AvaiContUsed) == "yes")
  if(arguments$ctrl_list_path == "NA" && !fam_file){
    logr::log_print("No control path defined and sample list is not a fam file")
    controls <- controls_1 %>%
      filter(tolower(BroadPhenotype) %in% tolower(BroadPhenotype_include))
    if(arguments$DetailedPhenotype_excl_regex != "NA"){
      logr::log_print("Applying regular expressoin to exclude specific detailed phenotypes")
      controls <- controls %>% filter(!str_detect(DetailedPhenotype, regex(arguments$DetailedPhenotype_excl_regex, ignore_case = TRUE))) 
    }
    logr::log_print(str(controls))
  } else if (fam_file) { #using a fam file
    logr::log_print(ctrl_list)
    controls <- data %>% 
      filter(tolower(sample_internal_name) %in% tolower(ctrl_list$V1))
    cases <- data %>% 
      filter(tolower(sample_internal_name) %in% tolower(case_list$V1))
  } else {#specific ctrl list pre-specified
    logr::log_print(ctrl_list <- fread(here("Input", arguments$ctrl_list_path), header = FALSE))
    controls <- controls_1 %>% 
      filter(tolower(sample_internal_name) %in% tolower(ctrl_list$V1))
  }
  logr::log_print(controls %>% select(BroadPhenotype) %>% table())
  
  # controls <- filtData %>% 
  #   filter(tolower(sample_internal_name) %!in% tolower(case_list$V1), seqGender %in% c("M", "F")) %>%
  #   filter(AvaiContUsed %in% c("yes", "Yes")) %>% 
  #   filter(tolower(BroadPhenotype) %in% tolower(BroadPhenotype_include)) %>%
  #   filter(!(tolower(DetailedPhenotype) %like any% c("%epilep%","%seizure%"))) 
  
  
  logr::log_print(paste0("Number of high quality controls in filtered dragen db: ", nrow(controls %>% distinct(sample_internal_name))))
  
  # The following code was gundi's original. More clever about finding relationships. Should be incorporated if possible.
  # controls <- filtData %>% 
  #     filter(tolower(sample_internal_name) %!in% tolower(case_list$V1), seqGender %in% c("M", "F")) %>%
  #     filter(AvaiContUsed %in% c("yes", "Yes")) %>% 
  #     filter(BroadPhenotype %in% c("healthy family member", "control", "kidney and urological disease",
  #                                  "amyotrophic lateral sclerosis", "obsessive compulsive disorder", 
  #                                  "dementia", "epilepsy", "liver disease", "other neurological disease",
  #                                  "other neuropsychiatric disease", "other neurodevelopmental disease", 
  #                                  "intellectual disability", "control mild neuropsychiatric disease") | 
  #                (BroadPhenotype %in% c("", "other") & str_detect(SubProject, "ALS"))) %>%  
  #     #            filter(BroadPhenotype %in% c("healthy family member", "control"))
  #     filter(!str_detect(DetailedPhenotype, regex("pulmo|fibrosis|lung", ignore_case = TRUE))) 
  
  if (!fam_file){
    dup2u <- controls %>% filter(duplicated(paste(sample_internal_name, TotalReads, sep="_"))| 
                                   duplicated(paste(sample_internal_name, TotalReads, sep="_"), fromLast=TRUE)) %>% 
      distinct(sample_internal_name, TotalReads, .keep_all = TRUE)
  
    dup3u <- controls %>% filter((duplicated(sample_internal_name) | duplicated(sample_internal_name, fromLast=TRUE)) & 
                                   !sample_internal_name %in% dup2u$sample_internal_name) %>% 
      group_by(sample_internal_name) %>%
      top_n(1, CCDSBasesCov10X) %>% 
      distinct(sample_internal_name, .keep_all = TRUE)
  
    controlsU <- controls %>% filter(!(duplicated(sample_internal_name)| 
                                         duplicated(sample_internal_name, fromLast=TRUE)))
  
    controlsT <- bind_rows(controlsU, dup2u, dup3u)
    controls <- controlsT
  }
  # cases in seqdb excluded after filtering 
  if(arguments$ctrl_list_path != "NA"){
    excluded_ctrl_seq_data <- data %>% filter(tolower(sample_internal_name) %in% tolower(ctrl_list$V1),
                                            tolower(sample_internal_name) %!in% tolower(controls$sample_internal_name))
    logr::log_print(paste0("Number of ctrls excluded by filtering: ", nrow(excluded_ctrl_seq_data %>% distinct(sample_internal_name))))#not sure if distinct is the best command here
    if (nrow(excluded_ctrl_seq_data) > 0){
      excluded_ctrl_seq_data$DetailedPhenotype <- NULL
      excluded_ctrl_seq_data$experimentID <- NULL
      excluded_ctrl_seq_data$sample_external_name <- NULL
      excluded_ctrl_seq_data$sample_aka <- NULL
      excluded_ctrl_seq_data$subproject_id <- NULL
      excluded_ctrl_seq_data$SubProject <- NULL
      excluded_ctrl_seq_data$status_date <- NULL    
      excluded_ctrl_seq_data$submission_date <- NULL
      excluded_ctrl_seq_data$AgeAtCollection <- NULL
      excluded_ctrl_seq_data$SelfDeclEthnic <- NULL
      excluded_ctrl_seq_data$SelfDeclEthnicDetail <- NULL
      write.csv(excluded_ctrl_seq_data, 
                file = here(paste0("Data/", time_case_prefix,"excluded_ctrls.csv")),row.names=FALSE)
      exclusion_reason <- excluded_ctrl_seq_data %>%
      mutate("Problem" = case_when(
        !(sample_status == "In DragenDB") ~ "sample_status",
        !(sample_type %in% c("Exome", "Genome")) ~ "sample_types",
        (tolower(sample_internal_name) %in% tolower(exclude_list)) ~ "exclude_list",
        !(ContaminationPercentage <= ContaminationPercent_var | (ContaminationPercentage <= 15 & MeanCoverage > 200)) ~ "ContaminationPercentage",
        !(CCDSBasesCov10X >= CCDSBasesCov10X_var) ~ "CCDSBasesCov10X",
        !(tolower(capture_kit) %in% tolower(acceptable_capture_kits)) ~ "capture_kit",
        !(DBSNPOverlapSNVs >= DBSNPOverlapSNVs_var) ~ "DBSNPOverlapSNVs",
        !(DBSNPOverlapIndels >= DBSNPOverlapIndels_var) ~ "DBSNPOverlapIndels",
        (sample_internal_name %in% gnomAD$V1) ~ "gnomAD",
        !(toupper(seqGender) %in% c("M", "F")) ~ "gender",
        !(tolower(sample_internal_name) %!in% tolower(case_list$V1)) ~ "in case list",
        !(tolower(AvaiContUsed) == "yes") ~ "not available for continued use")) %>%
        select(sample_internal_name, Problem)
      write.csv(exclusion_reason,
                file = here(paste0("Data/", time_case_prefix,"reason for_excluded_ctrls.csv")),row.names=FALSE)
    }  
  } else if (fam_file) { #controls not in dragen
    controls_not_in_dragen <- ctrl_list %>% filter(tolower(V1) %!in% tolower(data$sample_internal_name))
    if(nrow(controls_not_in_dragen)){
    logr::log_print(paste0("There are controls absent from dragen: ", controls_not_in_dragen$V1))
    write.table(file = here(paste0("Data/", time_case_prefix,"controls_absent_from_dragen.txt")), 
                x = controls_not_in_dragen$V1, row.names = FALSE, col.names = FALSE, quote = FALSE) 
    }
  }
  
  #########################################################################
  
  # Option 1:
  ######### change this ############
  # caseData <- fread(paste0(sample_path, "sample_result_16940084_1561752175.csv")) # case only csv file downloaded from https://sequence.igm.cumc.columbia.edu/
  # 
  # if (!is.null(caseData$capture_kit)) {
  #     caseData <- caseData %>% 
  #         mutate(capture_kit =  case_when(capture_kit == ""|is.na(capture_kit) ~ exomeKit, TRUE ~ capture_kit))
  #     
  # } else {
  #     caseData <- caseData %>% 
  #         mutate(capture_kit =  case_when(sample_type == "Genome" ~ "Roche", TRUE ~ exomeKit))
  #     
  # }
  # 
  # # make sure to use same filters used for controls above
  # filtCaseData <- caseData %>% 
  #     filter(sample_status == "In DragenDB") %>% 
  #     filter(sample_type %in% c("Exome", "Genome")) %>%
  #     filter(ContaminationPercentage <= 2 | (ContaminationPercentage <= 15 & MeanCoverage > 200)) %>% 
  #     filter(CCDSBasesCov10X >= 90) %>% # note: new threshold
  #     filter(exomeKit %in% c("65MB", "Roche", "RocheV2", "IDTERPv1", "IDTERPv1Plus", "IDTxGEN", 
  #                            "IDTERPv1mtDNA", "AgilentCRE", "AgilentV4", "Agilentv5", "AgilentV5", 
  #                            "AgilentV5UTR", "AgilentV6", "MedExome", "N/A")) %>% 
  #     filter(DBSNPOverlapSNVs >= 0.85 & DBSNPOverlapIndels >= 0.8) %>% # new
  #     filter(!sample_internal_name %in% gnomAD$CHGVID) 
  # 
  # # or Option 2:
  # filtCaseData <- filtData %>% filter(SubProject %in% c("IPF")) 
  
  if (!fam_file){
    filtCaseData <- filtData %>% filter(tolower(sample_internal_name) %in% tolower(case_list$V1))
    logr::log_print(paste0("Number of cases in filtered dragen db: ", nrow(filtCaseData)))
  
    # remove duplicates
    dup2u <- filtCaseData %>% filter(duplicated(paste(sample_internal_name, TotalReads, sep="_"))| duplicated(paste(sample_internal_name, TotalReads, sep="_"), fromLast=TRUE)) %>% distinct(sample_internal_name, TotalReads, .keep_all = TRUE)
    logr::log_print(paste0("Number of cases in filtered dragen db with duplicates 2: ", nrow(dup2u)))
  
    dup3u <- filtCaseData %>% filter((duplicated(sample_internal_name) | duplicated(sample_internal_name, fromLast=TRUE)) & !sample_internal_name %in% dup2u$sample_internal_name) %>% 
      group_by(sample_internal_name) %>% 
      top_n(1, CCDSBasesCov10X) %>% 
      distinct(sample_internal_name, .keep_all = TRUE)
  
    logr::log_print(paste0("Number of cases in filtered dragen db with duplicates 3: ", nrow(dup3u)))
                         
    filtCaseDataU <- filtCaseData %>% filter(!(duplicated(sample_internal_name)| duplicated(sample_internal_name, fromLast=TRUE)))
    logr::log_print(paste0("Number of cases in filtered dragen db post duplicate removal: ", nrow(filtCaseDataU)))
  
    filtCaseDataT <- bind_rows(filtCaseDataU, dup2u, dup3u)
    filtCaseData <- filtCaseDataT
    logr::log_print(paste0("Number of cases in filtered dragen db pre-sex mismatch removal: ",nrow(filtCaseData)))
  
    # only include samples without sex mismatch
    filtCaseData2 <- filtCaseData %>%
      filter(SelfDeclGender == "Unknown" | 
                 (SelfDeclGender == seqGender) |
                 (SelfDeclGender == "M" &
                      seqGender == "Ambiguous" & XYAvgCovRatio < 20))
  
    logr::log_print(paste0("Number of cases in filtered dragen db post-sex mismatch removal: ",nrow(filtCaseData2)))
  
    # cases in seqdb excluded after filtering 
    excluded_case_seq_data <- data %>% filter(tolower(sample_internal_name) %in% tolower(case_list$V1),
      tolower(sample_internal_name) %!in% tolower(filtCaseData2$sample_internal_name))
    logr::log_print(paste0("Number of cases excluded by filtering: ", nrow(excluded_case_seq_data %>% distinct(sample_internal_name))))
    cases <- filtCaseData2
  
    if (nrow(excluded_case_seq_data) > 0){
    excluded_case_seq_data$DetailedPhenotype <- NULL
    excluded_case_seq_data$experimentID <- NULL
    excluded_case_seq_data$sample_external_name <- NULL
    excluded_case_seq_data$sample_aka <- NULL
    excluded_case_seq_data$subproject_id <- NULL
    excluded_case_seq_data$SubProject <- NULL
    excluded_case_seq_data$status_date <- NULL    
    excluded_case_seq_data$submission_date <- NULL
    excluded_case_seq_data$AgeAtCollection <- NULL
    excluded_case_seq_data$SelfDeclEthnic <- NULL
    excluded_case_seq_data$SelfDeclEthnicDetail <- NULL
    write.csv(excluded_case_seq_data, 
              file = here(paste0("Data/", time_case_prefix,"excluded_cases.csv")),row.names=FALSE)
    exclusion_reason <- excluded_case_seq_data %>%
    mutate("Problem" = case_when(
      !(sample_status == "In DragenDB") ~ "sample_status",
      !(sample_type %in% c("Exome", "Genome")) ~ "sample_types",
      (tolower(sample_internal_name) %in% tolower(exclude_list)) ~ "exclude_list",
      !(ContaminationPercentage <= ContaminationPercent_var | (ContaminationPercentage <= 15 & MeanCoverage > 200)) ~ "ContaminationPercentage",
      !(CCDSBasesCov10X >= CCDSBasesCov10X_var) ~ "CCDSBasesCov10X",
      !(tolower(capture_kit) %in% tolower(acceptable_capture_kits)) ~ "capture_kit",
      !(DBSNPOverlapSNVs >= DBSNPOverlapSNVs_var) ~ "DBSNPOverlapSNVs",
      !(DBSNPOverlapIndels >= DBSNPOverlapIndels_var) ~ "DBSNPOverlapIndels",
      (sample_internal_name %in% gnomAD$V1) ~ "gnomAD",
      !(SelfDeclGender == "Unknown" | (SelfDeclGender == seqGender)) ~ "SelfDeclGender/seqGender",
      !(SelfDeclGender == "M" & seqGender == "Ambiguous" & XYAvgCovRatio < 20) ~ "seqGender/XYAvgCovRatio" )) %>%
      select(sample_internal_name, Problem)
    write.csv(exclusion_reason, 
              file = here(paste0("Data/", time_case_prefix,"reason_for_excluded_cases.csv")),row.names=FALSE)
  }
  }
  # or Option 3:
  
  
  # 1) Family ID: specify a family ID or use the same value as Individual ID to indicate this sample is being used as a non family sample
  # 2) Individual ID: sample ID
  # 3) Paternal ID: 0
  # 4) Maternal ID: 0
  # 5) Sex: 1=male, 2=female
  # 6) Phenotype: 1=control, 2=case
  # 7) Sample Type
  # 8) Capture Kit
  
  
  samples <- rbind(cases, controls)
  if (!fam_file) {
    samples <- samples[!duplicated(samples$sample_internal_name),]
    samples <- samples %>% filter(tolower(sample_internal_name) %!in% tolower(exclude_list)) # there is something wrong with these 2 samples
    logr::log_print(paste0("Number of case/ctrl samples after excluding duplicates and excluding absolute exclusion list : ",nrow(samples)))
  }
  samplesF <- samples %>% mutate(FamID = sample_internal_name, IndivID = sample_internal_name,
                                 PatID = 0, MatID = 0,
                                 Sex = case_when(seqGender == "F" ~ 2, TRUE ~ 1),
                                 Pheno = case_when(sample_internal_name %in% cases$sample_internal_name ~ 2, TRUE ~ 1),
                                 SampleType = case_when(sample_type == "Genome" ~ "Genome_As_Fake_Exome", TRUE ~ sample_type),
                                 CaptureKit = capture_kit) %>% select(FamID, IndivID, PatID, MatID, Sex, Pheno, SampleType, CaptureKit)
  write.table(samplesF, file = here(paste0("Data/", time_case_prefix, "CasesAndControlsAll.ped.txt")), sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  # *******
  # write seq files
  # ******
  controls$DetailedPhenotype <- NULL
  controls$experimentID <- NULL
  controls$sample_external_name <- NULL
  controls$sample_aka <- NULL
  controls$subproject_id <- NULL
  controls$SubProject <- NULL
  controls$status_date <- NULL    
  controls$submission_date <- NULL
  controls$AgeAtCollection <- NULL
  controls$SelfDeclEthnic <- NULL
  controls$SelfDeclEthnicDetail <- NULL
  write.csv(controls, 
            file = here(paste0("Data/", time_case_prefix,"ctrl.seq.csv")),row.names=FALSE)
  logr::log_print(paste0("Number of cases samples included : ",nrow(samplesF %>% filter(Pheno == 2))))
  
  cases$DetailedPhenotype <- NULL
  cases$experimentID <- NULL
  cases$sample_external_name <- NULL
  cases$sample_aka <- NULL
  cases$subproject_id <- NULL
  cases$SubProject <- NULL
  cases$status_date <- NULL    
  cases$submission_date <- NULL
  cases$AgeAtCollection <- NULL
  cases$SelfDeclEthnic <- NULL
  cases$SelfDeclEthnicDetail <- NULL
  write.csv(cases, 
            file = here(paste0("Data/", time_case_prefix,"case.seq.csv")),row.names=FALSE)
  logr::log_print(paste0("Number of ctrl samples included : ",nrow(samplesF %>% filter(Pheno == 1))))
  
  write.csv(rbind(cases,controls), file = here(paste0("Data/", time_case_prefix,"full_data.seq.csv")),row.names=FALSE)
  logr::log_print(paste0("Exclude ",nrow(excluded_case_seq_data), " cases which are in DragenDB due to filtering"))
}, error = function(e) {
  message("Caught an error in remainder of code: ", e$message)
})

tryCatch({
  logr::log_print("Finished")
  logr::log_close()
}, error = function(e) {
  # Handle the error here
  print(paste("An error occurred while closing log file:", e$message))
})
