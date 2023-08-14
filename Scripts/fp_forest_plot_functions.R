library(tidyverse)
library(data.table)
library(DescTools)

#'returns combined dataframes with sample and genotype data for a cluster
cluster_sample_and_genotype_list_fp<-function(git_root,resolution_var, min_case, model, nclust_overide=NULL,matrix_geno = "geno"){

  print(sprintf("Loading in %s project, %s model, %i min_case, %s resolution", git_root, model, min_case, resolution_var))#, nclust_overide=nclust_overide, 
  
  sample_df <- cluster_sample_list_fp(git_root,resolution_var, min_case, model, nclust_overide=nclust_overide)
  nclust <- unique(sample_df$cluster)
  
  print("In cluster_sample_and_genotype_list, done sample df")#, nclust_overide=nclust_overide, 
  if(matrix_geno=="geno"){
    genotype_df<-cluster_genotype_df_fp(git_root,resolution_var, model, nclust)
  } else if(matrix_geno=="matrix"){
    genotype_df<-cluster_matrix_df_fp(git_root,resolution_var, model, nclust)
  } 
  genotype_df$model <- model
  
  
  return_list<-list()
  return_list$genotype_df <- genotype_df
  return_list$sample_df <- sample_df
  return_list$nclust<- nclust
  return(return_list)
  
}

#'returns data frame of sample files
cluster_genotype_df_fp<-function(git_root,resolution_var, model, nclust){
  
  # *********
  # List of genotype file locations
  # *********
  genotype_location_vector<-cluster_genotype_location_vector_fp(git_root,resolution_var, model, nclust,matrix_geno = "geno")
  genotype_list <- lapply(genotype_location_vector,return_genotype_with_cluster_fp,model)
  # temp<-return_genotype_with_cluster(genotype_location_vector[1], model,epi_phenotype_var = epi_phenotype_var )#TODO return and fix this
  genotype_df <- data.frame()
  # *********
  # loop and combine
  # *********
  # genotype_df <- do.call(rbind, genotype_list)
  for (i in 1:length(genotype_list)){
    # print(i)
    temp_2_df<-as.data.frame(genotype_list[i])
    temp_2_df$cluster<-nclust[i]
    genotype_df <- rbind(genotype_df,temp_2_df)
    # genotype_df$cluster <- nclust[i]
  }
  
  return(genotype_df)
} 

#' returns combined genotype file for clusters
cluster_matrix_df_fp<-function(git_root,resolution_var, model, nclust){
  
  # *********
  # List of genotype file locations
  # *********
  genotype_location_vector<-cluster_genotype_location_vector_fp(git_root,resolution_var, model, nclust,matrix_geno = "matrix")
  genotype_list <- lapply(genotype_location_vector,return_genotype_with_cluster_fp,model, matrix_geno = "matrix")
  # temp<-return_genotype_with_cluster(genotype_location_vector[1], model,epi_phenotype_var = epi_phenotype_var )#TODO return and fix this
  genotype_df <- data.frame()
  # *********
  # loop and combine
  # *********
  # genotype_df <- do.call(rbind, genotype_list)
  for (i in 1:length(genotype_list)){
    # print(i)
    temp_2_df<-as.data.frame(genotype_list[i])
    if(nrow(temp_2_df)==0) next
    temp_2_df$cluster<-nclust[i]
    genotype_df <- rbind(genotype_df,temp_2_df)
    # genotype_df$cluster <- nclust[i]
  }
  
  return(genotype_df)
} 

#' returns list of genotype file paths
cluster_genotype_location_vector_fp<-function(git_root,resolution_var, model, nclust,matrix_geno = "geno"){
  genotype_location_vector <-character(length=length(nclust))
  if(matrix_geno=="geno"){
    for(i in 1:length(nclust)){
      # print(i)
      if(!(model %like% "%ecessiv%")){
        print(genotype_location_vector[i] <- list.files(paste0(git_root,"Results/Collapsing/LClust_res_", resolution_var, "_cluster_", nclust[i], "_FlashColl_07/", model), "_genotypes.csv$", full.names = TRUE))
      } else {
        print(genotype_location_vector[i] <-  list.files(paste0(git_root,"Results/Collapsing/LClust_res_", resolution_var, "_cluster_", nclust[i], "_FlashColl_07/", model), "_comphet.csv$", full.names = TRUE))
      }
    }
  } else if (matrix_geno=="matrix"){
    for(i in 1:length(nclust)){
      print(genotype_location_vector[i] <- list.files(paste0(git_root,"Results/Collapsing/LClust_res_", resolution_var, "_cluster_", nclust[i], "_FlashColl_07/", model), "matrix.txt.gz$", full.names = TRUE))
    }
  }
  return(genotype_location_vector)
}

return_genotype_with_cluster_fp<-function(genotype_path,model,matrix_geno = "geno"){
  # print(paste("Loading genotype file... ", genotype_path, sep = ""))
  # if(matrix_geno=="geno"){
  temp_df<-atav_load_genotype_fp(genotype_path,model,matrix_geno=matrix_geno)
  # } else if (matrix_geno == "matrix"){
  #   temp_df<-atav_load_matrix_fp(genotype_path,model)
  # }
  print(paste("Done loading genotype file... ", genotype_path, sep = ""))
  temp_df$file_location <- rep(genotype_path,nrow(temp_df))
  return(temp_df)
}

#' Load genotype file and update it for code. THIS CAN BE DELETED IF NO ERRORS IN CURRENT MATRIX LOADING CODE.
# atav_load_matrix <- function(atav_path){
#   data <- read_delim(atav_path, delim = "\t")
#   
#   x<-function(i,data_matrix){
#     sample_name <- names(data_matrix)[i]
#     genes_positive <- data_matrix[data_matrix[,i]==1,1]
#     if(nrow(genes_positive)>0) return(data.frame("Sample.Name" = rep(sample_name,nrow(genes_positive)), "Gene.Name" = paste0("'",genes_positive$`sample/gene`,"'")))
#     return(NULL)
#   }
#   data_frame_list <- lapply(2:length(names(data)),x,data)
#   atav_df <- do.call(rbind, data_frame_list)
#   
# }


#' Load genotype file and update it for code
atav_load_genotype_fp <- function(atav_path,model,matrix_geno = "geno"){
  # *********
  # read in atav data
  # *********
  # print(atav_path)
  print(paste("Loading in genetype file... ",atav_path, sep = ""))#
  # atav_df <- read.csv(file=atav_path, header=TRUE, sep=",")
  if(matrix_geno == "geno"){
    if(!(model %like% "%ecessiv%")){
      # atav_df <- fread(atav_path, showProgress=TRUE)
      atav_df <- (fread(atav_path, showProgress=TRUE))# %>% filter((Effect %in% c(FUNCTIONAL_EFFECTS,syn_effects)) ) 
      # atav_df <- atav_df[,names(atav_df) %like any% c("Variant%","%Effect%","%Transcript%","%Polyphen%","%gnomAD%","%Gene Name%","Gene Symbols","%Sample%","%ClinVar%","%pLI%","%MTR%","%LIMBR%","REVEL","%PEXT%","%HGVS%")]
    } else {
      # atav_df <- read.csv(atav_path) 
      atav_df <-  fread(atav_path, fill = TRUE, sep = ",",showProgress = TRUE)
    }
    atav_extra_rows <- atav_add_rows_for_multiple_transcripts_fp(atav_df)
    print("Finished adding extra rows in non-recessive model")
  } else if(matrix_geno=="matrix"){
    data <- read_delim(atav_path, delim = "\t")
    
    x<-function(i,data_matrix){
      sample_name <- names(data_matrix)[i]
      genes_positive <- data_matrix[data_matrix[,i]==1,1]
      if(nrow(genes_positive)>0) return(data.frame(Sample.Name = rep(sample_name,nrow(genes_positive)), Gene.Name = paste0("'",genes_positive$`sample/gene`,"'")))
      return(NULL)
    }
    if(nrow(data) > 0){
      data_frame_list <- lapply(2:length(names(data)),x,data)
      atav_df <- do.call(rbind, data_frame_list)
      atav_extra_rows <- atav_df
    } else {
      atav_extra_rows <- data.frame(matrix(ncol = 2, nrow = 0))
      names(atav_extra_rows) <- c("Sample.Name","Gene.Name")
    }
    
  }
  
  print("Now returning atav_df")
  return(atav_extra_rows)
}

#'returns data frame of sample files
cluster_sample_list_fp<-function(git_root,resolution_var, min_case, model, nclust_overide=NULL){
  if(!is.null(nclust_overide)) {
    sample_location_vector<-cluster_sample_location_vector_fp(git_root,resolution_var, min_case, model,nclust_overide = nclust_overide)
  } else {
    sample_location_vector<-cluster_sample_location_vector_fp(git_root,resolution_var, min_case,model)
  }
  
  nclust <- sample_location_vector$nclust
  sample_df<- cluster_sample_df_fp(sample_location_vector$sample_location_vector, nclust)
  return(sample_df)
}

#' returns sample location vector for cluster algorithm
cluster_sample_location_vector_fp <- function(git_root,resolution_var, min_case,model,nclust_overide = NA){
  logr::log_print(cl_path <- list.files(path = paste0(git_root,"Results/Plots"), pattern = paste0("res_", resolution_var,"_cluster_sizes.txt$"), full.names = TRUE))
  if(!is.na(nclust_overide)) {
    nclust <- nclust_overide
  } else if(min_case > 0) {
    cl_sizes <- fread(cl_path)
    logr::log_print(nclust <- (cl_sizes %>% filter(case >= min_case, control>=min_case))$cluster)
  } else if(min_case == 0) {
    cl_sizes <- fread(cl_path)
    nclust <- cl_sizes$cluster
  }
  sample_location_vector <-character(length=length(nclust))
  for(i in 1:length(nclust)){
    sample_location_vector[i]<- list.files(path = paste0(git_root,"Results/Collapsing/LClust_res_",resolution_var,"_cluster_",nclust[i],"_FlashColl_07/",model), 
                                           pattern = "*_existing.sample.txt$", full.names = TRUE)

  }
  a<-list()
  a$sample_location_vector <- sample_location_vector
  a$nclust <- nclust
  return(a)
}

#' returns sample location vector for cluster algorithm
cluster_sample_df_fp <-function(sample_location_vector, nclust){
  x <- function(sample_path){
    # print(paste("Loading sample file... ", sample_path, sep = ""))
    temp_df<-atav_load_sample_fp(sample_path)
    temp_df$file_location <- rep(sample_path,nrow(temp_df))
    return(temp_df)
  }
  sample_list <-lapply(sample_location_vector,x)
  sample_df <- data.frame()
  for (i in 1:length(sample_list)){
    temp_sample_df <-as.data.frame(sample_list[i])
    temp_sample_df$cluster <- nclust[i]
    sample_df <- rbind(sample_df,temp_sample_df)
  }
  return(sample_df)
}

#' Load sample file and update it for code
atav_load_sample_fp <- function(sample_path){
  # ****
  # read in ped file
  # ****
  print(paste("Reading sample file... ", sample_path, sep = ""))
  # sample_df <- read.table(file=sample_path, header = FALSE, sep = "\t")
  sample_df <- fread(sample_path)
  
  return(sample_df)
}  

#' code here potentially for splitting genotypes entries into multiple rows based on 
atav_add_rows_for_multiple_transcripts_fp <- function(genotype_df){
  df_names <- names(genotype_df)
  df_names_updated <- gsub(" ","\\.",df_names)
  names(genotype_df) <- df_names_updated
  
  # print("In atav_add_rows_for_multiple_transcripts 1")
  # *********
  # split into multi and single transcript dfs
  # *********
  genotype_multi_transcript <- (genotype_df %>% filter(`Consequence.annotations:.Effect|Gene|Transcript|HGVS_c|HGVS_p|Polyphen_Humdiv|Polyphen_Humvar` %like any% c("%,%")))  %>% arrange(desc(Variant.ID))
  # genotype_multi_transcript$Transcript <- genotype_multi_transcript$Transcript.Stable.Id
  
  genotype_single_transcript <- (genotype_df %>% filter(!(`Consequence.annotations:.Effect|Gene|Transcript|HGVS_c|HGVS_p|Polyphen_Humdiv|Polyphen_Humvar` %like any% c("%,%")))) %>% arrange(desc(Variant.ID)) 
  # print("In atav_add_rows_for_multiple_transcripts 2")
  # genotype_single_transcript$Transcript <- genotype_single_transcript$Transcript.Stable.Id
  # genotype_single_transcript$`Consequence.annotations:.Effect|Gene|Transcript|HGVS_c|HGVS_p|Polyphen_Humdiv|Polyphen_Humvar` <- NULL
  
  # *********
  # take multitranscript and separate into multiple rows
  # *********
  # head(genotype_multi_split_df)
  genotype_multi_split_df <- separate_rows(genotype_multi_transcript,Variant.ID,Sample.Name, Gene.Name, UpToDate.Gene.Name, `Consequence.annotations:.Effect|Gene|Transcript|HGVS_c|HGVS_p|Polyphen_Humdiv|Polyphen_Humvar`, sep=",") %>% arrange(desc(Variant.ID))
  
  # print("In atav_add_rows_for_multiple_transcripts 3")
  # *********
  # update effect, gene.name, transcript, code change, and polyphen to reflect this transcript change
  # *********
  genotype_multi_split_update_columns_df <- genotype_multi_split_df %>% separate(`Consequence.annotations:.Effect|Gene|Transcript|HGVS_c|HGVS_p|Polyphen_Humdiv|Polyphen_Humvar`, c("Effect", "Gene.Name", "Transcript.Stable.Id", "HGVS_c", "HGVS_p", "Polyphen.Humdiv.Score", "Polyphen.Humvar.Score"), sep = "\\|", remove = FALSE)
  # this function worked prior to 1/2022. There was an ATAV change where the (CCDS) was removed from the column name
  # genotype_multi_split_update_columns_df <- genotype_multi_split_df %>% separate(`Consequence.annotations:.Effect|Gene|Transcript|HGVS_c|HGVS_p|Polyphen_Humdiv|Polyphen_Humvar`, c("Effect", "Gene.Name", "Transcript.Stable.Id", "HGVS_c", "HGVS_p", "Polyphen.Humdiv.Score.(CCDS)", "Polyphen.Humvar.Score.(CCDS)"), sep = "\\|", remove = FALSE)
  # print("In atav_add_rows_for_multiple_transcripts 4")
  genotype_multi_split_update_columns_df$Gene.Name <- paste0("'",genotype_multi_split_update_columns_df$Gene.Name,"'")
  # print("In atav_add_rows_for_multiple_transcripts 5")
  
  # *********
  # combine but keep only one change per position, gene and sample. Added in multi and single first to ensure distinct would keep most damaging variant in case of a tie
  # *********
  # genotype_multi_transcript$`Consequence.annotations:.Effect|Gene|Transcript|HGVS_c|HGVS_p|Polyphen_Humdiv|Polyphen_Humvar` <- NULL
  
  updated_genotype_df_temp <- rbind(genotype_multi_transcript,rbind(genotype_single_transcript,genotype_multi_split_update_columns_df)) 
  updated_genotype_df <- updated_genotype_df_temp[!duplicated(updated_genotype_df_temp[ , c("Variant.ID","Gene.Name", "Sample.Name")]), ]
  # Debugging
  # setdiff(names(genotype_single_transcript),names(genotype_multi_split_update_columns_df))
  # setdiff(names(genotype_multi_split_update_columns_df),names(genotype_single_transcript))# in genotype_multi_split_update_columns_df but not in genotype_single_transcript
  return(updated_genotype_df) #i think I don't need to worry about various gene name issues.
  genotype_original <- updated_genotype_df
  #if gene name doesn't match uptdodate gene name
  non_matches <- updated_genotype_df$Gene.Name != updated_genotype_df$UpToDate.Gene.Name
  #if gene name in table, use corresponding uptodate gene name
  x <- function(gene_name, gene_map_file){
    if(gene_name %in% gene_map_file$dragendb_gene){
      return(gene_map_file$uptodate_gene[gene_map_file$dragendb_gene == gene_name])
    } else {
      return(gene_name)
    }
  }
  updated_genotype_df$UpToDate.Gene.Name[non_matches] <- sapply(updated_genotype_df$Gene.Name[non_matches],x,gene_map_file)
  return(updated_genotype_df)
  # "FAM46D" %in% gene_map_file$dragendb_gene
  # gene_map_file$uptodate_gene[gene_map_file$dragendb_gene == "FAM46D"]
  #if not, uptodate gene name should equal gene name
  new_names <- names(updated_genotype_df)
  new_names[!(new_names %in% df_names_updated)]
}

#' The following code generates Figure: forest plot showing odds ratio with fisher test with atav output filtered by gene group. The rows will be 1st. 
#' @param df_list 
#' @param phenotypes_to_compare_1 
#' @param phenotypes_to_compare_2
#' @param gene_list names must be defined or will error
#' @param condition_vector names must be defined or will error
#' @param effect_list names must be defined or will error
#' @param title_str title of analysis
#' @param mc_meth sets multiple comparison method. Default is 'fdr'
#' @param cmh_flag sets CMH flag. Default is TRUE
#' @param ... loo_cluster_af (loo_cluster_af not defined --> rep(NA, length(phenotypes_to_compare_1))), 
#'       pub_names (pub_names not defined --> names(condition vector) ), 
#'       maf_min_nonneuro_vector <- rep(c(NA), length(phenotypes_to_compare_1))
#'       maf_max_nonneuro_vector <- rep(c(NA), length(phenotypes_to_compare_1))
#'       maf_min_vector <- rep(c(NA), length(phenotypes_to_compare_1))
#'       maf_max_vector <- rep(c(NA), length(phenotypes_to_compare_1))
#' @return Returns
#' @examples
#' figure_or_forest(df_list, phenotypes_to_compare_1, phenotypes_to_compare_2,gene_list,condition_vector,effect_list,title_str, ..., mc_meth = "fdr", cmh_flag = TRUE)#need to fill this in
figure_or_forest_fp<-function(atav_output_1, sample_df, phenotypes_to_compare_1, phenotypes_to_compare_2,
                           gene_list,condition_vector,title_str, ...,
                           mc_meth = "fdr", test_flag = "cmh", cov_var = NULL, suppress_output = TRUE,
                           cores = 1, save_filtered_df_folder = NULL){
  opt_args <- list(...)
  # Required variables
  if(is.null(names(gene_list))){#sets LOO
    error("Need to name gene list")
  }
  if(is.null(names(condition_vector))){#makes sure conditions have names
    error("Need to name condition vector")
  }
  # unlist variables
  if(!is.null(opt_args$effect_list)){#sets LOO
    effect_list=opt_args$effect_list
  } else {
    if(suppress_output == FALSE) print("effect list not defined so setting to NA")
    effect_list<-rep(NA, length(phenotypes_to_compare_1))
    names(effect_list) <- "NA"
  }
  if(!is.null(opt_args$loo_cluster_af)){#sets LOO
    loo_cluster_af=opt_args$loo_cluster_af
  } else {
    if(suppress_output == FALSE) print("loo_cluster_af not defined so setting to NA")
    loo_cluster_af<-rep(NA, length(phenotypes_to_compare_1))
  }
  if(is.null(opt_args$pub_names)){#sets LOO
    if(suppress_output == FALSE) print("pub_names not defined so setting to condition vector")
    pub_names <- names(condition_vector)
  } else {
    pub_names<-opt_args$pub_names
  }
  # print(pub_names)
  # *************
  # set non neuro MAFs
  # *************
  if(is.null(opt_args$maf_min_nonneuro_vector)){#sets population based MAF filter for non neuro
    if(suppress_output == FALSE) print("maf_min_nonneuro_vector not defined so setting to NA")
    maf_min_nonneuro_vector <- rep(c(NA), length(phenotypes_to_compare_1))
  } else{
    maf_min_nonneuro_vector <- opt_args$maf_min_nonneuro_vector
  }
  if(is.null(opt_args$maf_max_nonneuro_vector)){#sets population based MAF filter for non neuro
    if(suppress_output == FALSE) print("maf_max_nonneuro_vector not defined so setting to NA")
    maf_max_nonneuro_vector <- rep(c(NA), length(phenotypes_to_compare_1))
  } else{
    maf_max_nonneuro_vector <- opt_args$maf_max_nonneuro_vector
  }
  # *************
  # set population MAFs
  # *************
  if(is.null(opt_args$maf_min_vector)){#sets population based MAF filter for non neuro
    if(suppress_output == FALSE) print("maf_min_vector not defined so setting to NA")
    maf_min_vector <- rep(c(NA), length(phenotypes_to_compare_1))
  } else{
    maf_min_vector <- opt_args$maf_min_vector
  }
  if(is.null(opt_args$maf_max_vector)){#sets population based MAF filter for non neuro
    if(suppress_output == FALSE) print("maf_max_vector not defined so setting to NA")
    maf_max_vector <- rep(c(NA), length(phenotypes_to_compare_1))
  } else{
    maf_max_vector <- opt_args$maf_max_vector
  }
  
  
  
  # atav_output_1<-df_list$genotype_df
  # sample_df <- df_list$sample_df
  array_column_names <-c("Effect","Gene_group","Condition", "Group", "Case w/ QV","Case w/o QV", "Case % w/ QV", "Ctrl w/ QV","Ctrl w/o QV", "Ctrl % w/ QV", "OR", "[CI CI]","P_uncorrected","CI_2p5","CI_97p5","Filter","Group_1","Group_2","Publication Name","OR_original")
  array_for_plot <- matrix(data="",nrow=length(condition_vector),ncol=length(array_column_names),byrow = TRUE)
  colnames(array_for_plot)<-array_column_names
  temp_function <- function(temp_str) if(temp_str=="ctrl") return("Control") else if (temp_str=="case") return("Case") else return(temp_str)
  # i<-1
  # cycle for loop
  
  # atav_filtered_df_list <- mcmapply(return_atav_df_for_filtered,
  #                              rep(list(atav_output_1),length(condition_vector)),
  #                              gene_var = gene_list,
  #                              effect_var = effect_list,
  #                              condition_var=condition_vector,
  #                              loo_cluster_af=loo_cluster_af,
  #                              maf_max_nonneuro=maf_max_nonneuro_vector,
  #                              maf_min_nonneuro=maf_min_nonneuro_vector,
  #                              maf_max=maf_max_vector,
  #                              maf_min=maf_min_vector,
  #                              mc.cores = numCores, SIMPLIFY = FALSE)
  # 
  parallel_test <-function(dummy_var, atav_list, cov_var_list = NULL){
    # print(pub_names[dummy_var])
    atav_filtered_df<-return_atav_df_for_filtered_repo(atav_list, gene_var = gene_list[[dummy_var]], effect_var = effect_list[[dummy_var]],
                                                  condition_var=condition_vector[dummy_var],
                                                  loo_cluster_af=loo_cluster_af[dummy_var], maf_max_nonneuro=maf_max_nonneuro_vector[dummy_var],
                                                  maf_min_nonneuro=maf_min_nonneuro_vector[dummy_var], maf_max=maf_max_vector[dummy_var],
                                                  maf_min=maf_min_vector[dummy_var], suppress_output=suppress_output)
    if(!is.null(save_filtered_df_folder)) fwrite(atav_filtered_df, paste0(output_directory,dummy_var,"_",title_str,"_sanity_check.csv",sep = ""),na="NA")
    
    
    return(switch(test_flag,
                  "fet" = atav_fisher_fp(atav_filtered_df, sample_df, c(phenotypes_to_compare_1[dummy_var],phenotypes_to_compare_2[dummy_var]),cmh_test=FALSE),
                  "cmh" = atav_fisher_fp(atav_filtered_df, sample_df, c(phenotypes_to_compare_1[dummy_var],phenotypes_to_compare_2[dummy_var]),cmh_test=TRUE),
                  "log" = josh_implement_cov_regr(atav_filtered_df, sample_df, cov_var_list))  )
    
  }
  test_list_list <- mclapply(1:length(phenotypes_to_compare_1),
                             parallel_test,
                             atav_output_1, cov_var_list = cov_var, mc.cores = cores)
  # test_list_list <- parallel_test(1,  atav_output_1)
  
  array_for_plot[, "Effect"] <- names(effect_list)
  array_for_plot[, "Gene_group"] <- names(gene_list)
  array_for_plot[, "Condition"] <- condition_vector
  array_for_plot[, "Group_1"] <- phenotypes_to_compare_1
  array_for_plot[, "Group_2"] <- phenotypes_to_compare_2
  array_for_plot[, "Filter"] <- names(condition_vector)
  array_for_plot[, "Publication Name"] <- pub_names
  
  for (i in 1:length(condition_vector)){
    # print(i)
    # atav_filtered_df<-atav_filtered_df_list[[i]]
    # atav_filtered_df<-return_atav_df_for_filtered(atav_output_1, gene_var = gene_list[[i]], effect_var = effect_list[[i]],
    #                                               condition_var=condition_vector[i],
    #                                               loo_cluster_af=loo_cluster_af[i], maf_max_nonneuro=maf_max_nonneuro_vector[i],
    #                                               maf_min_nonneuro=maf_min_nonneuro_vector[i], maf_max=maf_max_vector[i],maf_min=maf_min_vector[i], suppress_output=TRUE)
    # if(nrow(df_to_write <- atav_filtered_df %>% filter(epi_phenotype == phenotypes_to_compare_1[i])) > 0)
    if(!is.null(opt_args$save_folder)){
      saveRDS(atav_filtered_df, file = paste0(opt_args$save_folder, title_str, "_included_row_",as.character(i),".rds"))
      saveRDS(setdiff(atav_output_1, atav_filtered_df[,1:ncol(atav_output_1)]), file = paste0(opt_args$save_folder, title_str, "_excluded_row_",as.character(i),".rds"))
    }
    test_list <- test_list_list[[i]]
    # test_list <- switch(test_flag,"fet" = atav_fisher(atav_filtered_df, sample_df, c(phenotypes_to_compare_1[i],phenotypes_to_compare_2[i]),cmh_test=FALSE),
    #        "cmh" = atav_fisher(atav_filtered_df, sample_df, c(phenotypes_to_compare_1[i],phenotypes_to_compare_2[i]),cmh_test=TRUE),
    #        "log" = josh_implement_cov_regr(atav_filtered_df, sample_df, cov_var))
    atav_test_var<-test_list[[1]]
    atav_test_var$conf.low <- atav_test_var$conf.int[1]
    atav_test_var$conf.high <- atav_test_var$conf.int[2]
    
    two_by_two_var <- test_list[[2]]
    array_for_plot[i, "Group"] <- sprintf("%s vs. %s",temp_function(phenotypes_to_compare_1[i]),temp_function(phenotypes_to_compare_2[i]))
    array_for_plot[i, "Case w/ QV"] <- two_by_two_var[1,1]
    array_for_plot[i, "Case w/o QV"] <- two_by_two_var[1,2]
    array_for_plot[i, "Case % w/ QV"] <- sprintf("%.1f%s",100*two_by_two_var[1,1] / (two_by_two_var[1,1] + two_by_two_var[1,2]),"%")
    array_for_plot[i, "Ctrl w/ QV"] <- two_by_two_var[2,1]
    array_for_plot[i, "Ctrl w/o QV"] <- two_by_two_var[2,2]
    array_for_plot[i, "Ctrl % w/ QV"] <- sprintf("%.1f%s",100*two_by_two_var[2,1] / (two_by_two_var[2,1] + two_by_two_var[2,2]),"%")
    array_for_plot[i, "OR"] <- sprintf("%.3f",as.numeric(atav_test_var$estimate))
    array_for_plot[i, "OR_original"] <- array_for_plot[i, "OR"]
    array_for_plot[i, "[CI CI]"] <- sprintf("[%.1f %.1f]", atav_test_var$conf.low, atav_test_var$conf.high)
    array_for_plot[i, "CI_2p5"] <- sprintf("%.3f", atav_test_var$conf.low)
    array_for_plot[i, "CI_97p5"] <- sprintf("%.3f", atav_test_var$conf.high)
    array_for_plot[i, "P_uncorrected"] <- atav_test_var$p.value
    # array_for_plot[i, "Mutation"] <- names(condition[i]
  }
  
  fisher_df <- as.data.frame(array_for_plot, stringsAsFactors = FALSE) %>% map_df(rev)
  # print(sprintf("Correcting for multiple comparisons with %s", mc_meth))
  fisher_df$P <- sprintf("%.1e", p.adjust(as.numeric(as.character(fisher_df$P_uncorrected)), method = mc_meth, n = nrow(fisher_df)))
  fisher_df$title_str <- title_str
  fisher_df$test_flag <- test_flag
  return(fisher_df)
}

# *********
# run fisher exact on variants in all genes in atav_df
# *********
atav_fisher_fp <- function(atav_df, peds_df, phenotypes_to_compare,cmh_test = FALSE){
  # if cmh_true, return CMH results instead of atav_fisher
  if(cmh_test==TRUE && (length(unique(peds_df$cluster)) > 1) ){
    two_by_two_array <- atav_create_2_by_2(atav_df, peds_df, phenotypes_to_compare, cmh_test=TRUE)
    two_by_two <- atav_create_2_by_2(atav_df, peds_df, phenotypes_to_compare)
    # print("in cmh")
    # print(two_by_two_array)
    return(list((mantelhaen.test(two_by_two_array, exact = TRUE)),two_by_two))
  }
  # ---------------
  # if cmh_true, return CMH results instead of atav_fisher
  # ---------------
  two_by_two <- atav_create_2_by_2(atav_df, peds_df, phenotypes_to_compare)
  # print(two_by_two)
  HE_number <- two_by_two["Healthy","Exposed"]
  HN_number <- two_by_two["Healthy","Not_Exposed"]
  DN_number <- two_by_two["Disease","Not_Exposed"]
  DE_number <- two_by_two["Disease","Exposed"]
  
  # ---------------
  # calculate fisher exact
  # ---------------
  fisher_val<-fisher.test(two_by_two, y = NULL, alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)
  # fisher_val<-fisher.test(two_by_two, y = NULL, or = 1, alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)
  # print(fisher_val)
  fisher_val$two_by_two <- two_by_two
  
  # print("leaving fisher")
  # library(tidyverse)
  return(a<-list((fisher_val),two_by_two))
}

# *********
# create 2 by 2 for odds ratio and fisher calculations
# *********
atav_create_2_by_2 <- function(atav_df, peds_df, phenotypes_to_compare, cmh_test = FALSE){
  # logr::log_print("updated")
  # *********
  # return 2 by 2 array for CMH testing
  # *********
  if(cmh_test==TRUE){
    cluster_vector <- unique(peds_df$cluster)
    array_2_by_2 <- array(dim=c(2,2,length(cluster_vector)))
    for(i in 1:length(cluster_vector)){
      temp_2_by_2 <- atav_create_2_by_2(atav_df %>% filter(cluster==cluster_vector[i]), peds_df %>% filter(cluster==cluster_vector[i]),phenotypes_to_compare)
      array_2_by_2[,,i] <- temp_2_by_2
    }
    return(array_2_by_2)
  }
  # print("entering atav_create_2_by_2")
  
  # *********
  # number of disease exposed (DE) and number of disease not exposed (DN)
  # *********
  if (phenotypes_to_compare[1] == "case"){
    case_names_total <- peds_df %>% filter(V6==2) %>% select(V2)
  } else {
    stop("This needs to be debugged")
    case_names_total <- peds_df %>% filter(epi_phenotype==phenotypes_to_compare[1]) %>% select(V2)
  }
  # if(!(model %like% "%ecessiv%")){
    case_names_exposed_temp <- atav_df %>% filter(as.character(Sample.Name) %in% as.character(case_names_total$V2))
    case_names_exposed <- case_names_exposed_temp[!duplicated(case_names_exposed_temp[ , c("Sample.Name")]), ] 
  # } else {
  #   stop("not tested")
  #   case_names_exposed <- unique(atav_df %>% filter(as.character(Sample.Name) %in% as.character(case_names_total$V2)) %>% select(Sample.Name))
  # }

  DE_number <- nrow(case_names_exposed)
  DN_number <- nrow(case_names_total) - nrow(case_names_exposed)
  if(DN_number < 0){
    logr::log_print("bad match between case genotype and ped file")
    error()
  }
  
  if(length(temp_case_names <- case_names_exposed$Sample.Name[!(case_names_exposed$Sample.Name %in% case_names_total$V2 )]) > 0){
    logr::log_print("bad match between case genotype and ped file")
    error()
  }
  # *********
  # number of healthy exposed (HE) and number of healthy not exposed (HN)
  # *********
  if (phenotypes_to_compare[2] == "ctrl"){
    ctrl_names_total <- peds_df %>% filter(V6==1) %>% select(V2)
  } else {
    stop("This needs to be debugged")
    ctrl_names_total <- peds_df %>% filter(epi_phenotype==phenotypes_to_compare[2]) %>% select(V2)
  }
  
  # ctrl_names_exposed <- unique(atav_df %>% filter(as.character(Sample.Name) %in% as.character(ctrl_names_total$V2)) %>% select(Sample.Name))
  # ctrl_names_exposed <- atav_df %>% filter(as.character(Sample.Name) %in% as.character(ctrl_names_total$V2)) %>% distinct(Sample.Name)
  ctrl_names_exposed_temp <- atav_df %>% filter(as.character(Sample.Name) %in% as.character(ctrl_names_total$V2))
  ctrl_names_exposed <- ctrl_names_exposed_temp[!duplicated(ctrl_names_exposed_temp[ , c("Sample.Name")]), ] 
  
  # ctrl_names_exposed <- atav_df %>% filter(epi_phenotype==phenotypes_to_compare[2]) %>% distinct(Sample.Name)#changed on 8/31
  HE_number <- nrow(ctrl_names_exposed)
  HN_number <- nrow(ctrl_names_total) - nrow(ctrl_names_exposed)
  if(HN_number < 0){
    logr::log_print("bad match between ctrl genotype and ped file")
    error()
  }
  if(length(temp_ctrl_names <- ctrl_names_exposed$Sample.Name[!(ctrl_names_exposed$Sample.Name %in% ctrl_names_total$V2 )]) > 0){
    logr::log_print("bad match between ctrl genotype and ped file")
    error()
  }
  # *********
  # 
  # *********
  two_by_two <- matrix(c(DE_number,DN_number, HE_number,HN_number),nrow=2,ncol=2, byrow = TRUE)
  dimnames(two_by_two) = list(c("Disease", "Healthy"), c("Exposed", "Not_Exposed"))
  
  return(two_by_two)
}


#' shuffles variant case control labels
#' by_cluster = TRUE means that shuffling happens within each cluster, FALSE (default) indicates shuffles case/control labels independent of cluster assignment
return_shuffled_atav_cluster_fp <-function(dummy_var, atav_temp_s, sample_temp_s){
  geno_list <- list()
  sample_list <- list()
  cluster_vector <- unique(sample_temp_s$cluster)
  # df_list <- lapply(cluster_vector, function(x) return_shuffled_atav(1,atav_temp %>% filter(cluster==x)))
  for(i in 1:length(cluster_vector)){
    temp_list <- return_shuffled_atav_2_fp(1,
                                        atav_temp_s %>% filter(cluster==cluster_vector[i]),
                                        sample_temp_s %>% filter(cluster==cluster_vector[i]))
    geno_list[[i]] <- temp_list$genotype_df
    sample_list[[i]] <- temp_list$sample_df
  }
  permutation_atav <- do.call(rbind, geno_list)
  permutation_sample <- do.call(rbind, sample_list)
  
  df_list <- list(permutation_atav, permutation_sample)
  names(df_list) <- c("genotype_df","sample_df")
  return(df_list)
}

return_shuffled_atav_2_fp <-function(dummy_var, atav_temp_2, sample_temp_2){
  shuffle_index<-sample(length(sample_temp_2$V6))
  
  permutation_sample <- sample_temp_2
  # permutation_sample$epi_phenotype <- sample_temp_2$epi_phenotype[shuffle_index]
  permutation_sample$V6 <- sample_temp_2$V6[shuffle_index]
  
  permutation_atav <- atav_temp_2#[1:5,]
  
  new_index <- match(permutation_atav$Sample.Name,permutation_sample$V1)
  # permutation_atav$epi_phenotype <- permutation_sample$epi_phenotype[new_index]
  permutation_atav$Sample.Phenotype <- permutation_atav$Sample.Phenotype[new_index]
  
  df_list <- list(permutation_atav, permutation_sample)
  names(df_list) <- c("genotype_df","sample_df")
  return(df_list)
  
}

# *********
# Filters atav 
# *********
return_atav_df_for_filtered_repo<-function(atav_df, maf_max = NA,maf_min = NA,gene_var = NA,effect_var = NA,condition_var = "Unfiltered",
                                      phenotype_vector = NA,suppress_output = TRUE, singleton = FALSE, loo_cluster_af = NA,
                                      maf_min_nonneuro = NA, maf_max_nonneuro = NA, only_common = FALSE){
  if( !is.na(gene_var[1])){
    old_atav_rows <- nrow(atav_df)
    atav_df <- atav_df %>% filter(Gene.Name %in% gene_var)
    if(suppress_output==FALSE){
      print("Filtering for Genes")
      print(sprintf("From %i to %i",old_atav_rows,nrow(atav_df)))
    }
  }
  if(!is.na(effect_var[1])){
    old_atav_rows <- nrow(atav_df)
    atav_df <- atav_df %>% filter(Effect %in% effect_var)
    if(suppress_output==FALSE){
      print(paste("Including only effect ", effect_var, sep=""))
      print(sprintf("From %i to %i",old_atav_rows,nrow(atav_df)))
    }
    
  }
  if(condition_var != "Unfiltered"){
    old_atav_rows <- nrow(atav_df)
    atav_df <- atav_df %>% filter(eval(parse(text=condition_var)))
    if(suppress_output==FALSE){
      print(paste("Including only condition ", condition_var, sep=""))
      print(sprintf("From %i to %i",old_atav_rows,nrow(atav_df)))
    }
  }
  if(!is.na(phenotype_vector)){
    old_atav_rows <- nrow(atav_df)
    new_atav_rows <- nrow(atav_df <- atav_df %>% filter(epi_phenotype %in% phenotype_vector))
    if(suppress_output==FALSE){
      print(paste("Including only epi phenotypes ", phenotype_vector, sep=""))
      print(sprintf("From %i to %i",old_atav_rows,new_atav_rows))
    }
    
  }
  if(singleton==TRUE){
    print("Singletons only")
    old_atav_rows <- nrow(atav_df)
    # error("dont use this")
    temp <- atav_df %>% group_by(Variant.ID) %>% summarize(count_var=n())
    include_variant_ids <- (temp %>% filter(count_var==1))$Variant.ID
    new_atav_rows <- nrow(atav_df <- atav_df %>% filter(Variant.ID %in% include_variant_ids))
    if(suppress_output==FALSE){
      print("Including only singlestons")
      print(sprintf("From %i to %i",old_atav_rows,new_atav_rows))
    }
  }
  
  if(!is.na(loo_cluster_af) ){
    old_atav_rows <- nrow(atav_df)
    temp <- atav_df %>% group_by(Variant.ID,cluster_loo=cluster) %>% summarize(count_var=(n()-1))
    cluster_size <- df_list$sample_df %>% group_by(cluster_num=cluster) %>% summarize(cluster_size_var=n())
    temp$cluster_size_var <-numeric(length=nrow(temp))
    for(cluster_loop in unique(temp$cluster_loo)){
      temp[temp$cluster_loo==cluster_loop, "cluster_size_var"] <- (cluster_size %>% filter(cluster_num==cluster_loop))$cluster_size_var
      
    }
    
    temp <- temp %>% mutate(af = count_var/(2*cluster_size_var), var_id_cluster = paste(Variant.ID,"-",as.character(cluster_loo),sep=""))
    
    include_variant_ids <- (temp %>% filter(af <= loo_cluster_af))$var_id_cluster
    # atav_df <- atav_df %>% mutate(var_id_cluster = paste(Variant.ID,"-",as.character(clus),sep=""))
    # atav_df <- atav_df %>% mutate(var_id_cluster = paste(Variant.ID,"-",as.character(cluster_loo),sep="")) %>% filter(var_id_cluster %in% include_variant_ids)
    
    atav_df <- atav_df %>% mutate(var_id_cluster = paste(Variant.ID,"-",as.character(cluster),sep="")) %>% filter(var_id_cluster %in% include_variant_ids)
    if(suppress_output==FALSE){
      print("LOO cluster")
      print(sprintf("From %i",old_atav_rows))
      print(sprintf(" to %i",nrow(atav_df)))
    }
  }
  # ***********
  # MAF non neuro 
  # ***********
  if(!is.na(maf_min_nonneuro)){
    old_atav_rows <- nrow(atav_df)
    suppressWarnings(atav_df <- atav_df %>% rowwise() %>% mutate(josh_gnomad_nonneuro_af=max(c(gnomAD.Exome.non_neuro_afr_AF,gnomAD.Exome.non_neuro_amr_AF,gnomAD.Exome.non_neuro_asj_AF, gnomAD.Exome.non_neuro_eas_AF,gnomAD.Exome.non_neuro_sas_AF,gnomAD.Exome.non_neuro_fin_AF,gnomAD.Exome.non_neuro_nfe_AF),na.rm = TRUE)))
    atav_df[atav_df$josh_gnomad_nonneuro_af==-Inf,"josh_gnomad_nonneuro_af"]<-0
    atav_df <- atav_df %>% filter(josh_gnomad_nonneuro_af > maf_min_nonneuro)
    
    if(suppress_output==FALSE){
      print("Non neuro MAF min")
      print(sprintf("From %i",old_atav_rows))
      print(sprintf(" to %i",nrow(atav_df)))
    }
  }
  if( !is.na(maf_max_nonneuro)){
    old_atav_rows <- nrow(atav_df)
    suppressWarnings(atav_df <- atav_df %>% rowwise() %>% mutate(josh_gnomad_nonneuro_af=max(c(gnomAD.Exome.non_neuro_afr_AF,gnomAD.Exome.non_neuro_amr_AF,gnomAD.Exome.non_neuro_asj_AF, gnomAD.Exome.non_neuro_eas_AF,gnomAD.Exome.non_neuro_sas_AF,gnomAD.Exome.non_neuro_fin_AF,gnomAD.Exome.non_neuro_nfe_AF),na.rm = TRUE)))
    atav_df[atav_df$josh_gnomad_nonneuro_af==-Inf,"josh_gnomad_nonneuro_af"]<-0
    atav_df <- atav_df %>% filter((josh_gnomad_nonneuro_af) <= maf_max_nonneuro)
    
    if(suppress_output==FALSE){
      print("Non neuro MAF max")
      print(sprintf("From %i",old_atav_rows))
      print(sprintf(" to %i",nrow(atav_df)))
    }
  }
  # ***********
  # MAF global
  # ***********
  if(!is.na(maf_min)){
    old_atav_rows <- nrow(atav_df)
    print(sprintf("From %i",nrow(atav_df)))
    suppressWarnings(atav_df <- atav_df %>% rowwise() %>% mutate(josh_gnomad_af=max(c(gnomAD.Exome.afr_AF,gnomAD.Exome.amr_AF,gnomAD.Exome.asj_AF, gnomAD.Exome.eas_AF,gnomAD.Exome.sas_AF,gnomAD.Exome.fin_AF,gnomAD.Exome.nfe_AF),na.rm = TRUE)))
    atav_df[atav_df$josh_gnomad_af==-Inf,"josh_gnomad_af"]<-0
    atav_df <- atav_df %>% filter(josh_gnomad_af > maf_min)
    
    if(suppress_output==FALSE){
      print("MAF min")
      print(sprintf("From %i",old_atav_rows))
      print(sprintf(" to %i",nrow(atav_df)))
    }
  }
  if(!is.na(maf_max)){
    old_atav_rows <- nrow(atav_df)
    suppressWarnings(atav_df <- atav_df %>% rowwise() %>% mutate(josh_gnomad_af=max(c(gnomAD.Exome.afr_AF,gnomAD.Exome.amr_AF,gnomAD.Exome.asj_AF, gnomAD.Exome.eas_AF,gnomAD.Exome.sas_AF,gnomAD.Exome.fin_AF,gnomAD.Exome.nfe_AF),na.rm = TRUE)))
    atav_df[atav_df$josh_gnomad_af==-Inf,"josh_gnomad_af"]<-0
    atav_df <- atav_df %>% filter((josh_gnomad_af) <= maf_max)
    
    if(suppress_output==FALSE){
      print("MAF max")
      print(sprintf("From %i",old_atav_rows))
      print(sprintf(" to %i",nrow(atav_df)))
    }
  }
  
  if(only_common==TRUE){
    old_atav_rows <- nrow(atav_df)
    common_genes <- return_common_genes(atav_df)
    atav_df<-atav_df[atav_df$Gene.Name %in% common_genes,]
    if(suppress_output==FALSE){
      print("Common genes only")
      print(sprintf("From %i",old_atav_rows))
      print(sprintf(" to %i",nrow(atav_df)))
    }
  }
  return(atav_df)
}
