# Created by Josh Motelow 1/26/2022
# modified implementation of unbiased threshold analysis conceived in Stanley, K. E. et al. Causal Genetic Variants in Stillbirth. N Engl J Med NEJMoa1908753 (2020) doi:10.1056/NEJMoa1908753. numCores and mac_pc are global variables that I haven't figured out how to handle yet.

library(logr)
library(parallel)
library(cowplot)

#' Measure optimal p value for missense intolerance threshold. 
#' @param genotype_df atav genotypes data frame
#' @param sample_df atav peds file 
#' @param analysis_label_var text string indicating label of analysis. Used to mark saved permutations and load them back in 
#' @param permutation_directory directory for permutations
#' @param output_str save directory for figures
#' @param nperm number of permutations to either create or to load in for display figure
#' @param save_stat determines the type of statistic saved as well as type of analysis performed. Options are "or" where just the OR is saved. "stat" means that whole comparison including or, p value, and test statistic are all saved.
#' @param test_flag determines type of test that you're running. Options are "cmh" "fet" "log"
#' @param decile_label_var i think location on plot of decile labels. 
#' @param max_y_disp maximum height of plot
#' @param y_break_marks location of labels
#' @return saves plots 
#' @examples
#' figure_gene_intolerance_threshold(genotype_df,sample_df, paste0("neuro_refractory_loeuf_UrPTV_min_case_",min_case,"_resolution_",resolution_var),permutation_directory, output_str=output_str,nperm=1000,save_stat = "stat", test_flag = "cmh")
#'  figure_gene_intolerance_threshold(genotype_df,sample_df,paste0("neuro_refractory_loeuf_UrPTV_min_case_",min_case,"_resolution_",resolution_var),permutation_directory,output_str=output_str,nperm=100,perform_permutations = FALSE,save_stat = "stat_p", test_flag = "cmh", decile_label_var = 3.5,max_y_disp = 6.1)
figure_gene_misz_intolerance_threshold<-function(genotype_df, sample_df, analysis_label_var,permutation_directory, output_str, nperm = 100000,
                                                  save_stat = "stat", test_flag=NULL, decile_label_var = 5.5, max_y_disp = 10, y_break_marks = c(0,2.5,5)){
  logr::log_open(paste0(output_str,"/", 
                        paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_misz_permutation_analysis_logfile_",analysis_label_var, "_"),"logfile.log"))
  
  if(is.null(test_flag)){
    logr::log_print("Exit: need to define test flag")
    stop()
  } 
  # ****************
  # initialize variables
  # ****************
  misz_seq <- unique(sort(genotype_df$mis_z), decreasing = TRUE)
  logr::log_print(sprintf("Looking at %i genes with %i unique Mis Z scores", length(unique(genotype_df$Gene.Name)), number_of_bins <- length(misz_seq)))
  gene_list <- mclapply(1:number_of_bins,
                        function(x) return(unique((genotype_df %>% filter(mis_z >= misz_seq[x]))$Gene.Name)),mc.cores = numCores )
  names(gene_list)<-rep(sprintf("Greater than or equal to %.3f",misz_seq),1)
  
  generate_unbiased_threshold_analysis(genotype_df, sample_df, analysis_label_var,permutation_directory, output_str, nperm, gene_list,misz_seq,
                                       save_stat, test_flag, decile_label_var, max_y_disp, y_break_marks, plot_type = "misz")
  
}

#' Measure optimal p value for loeuf intolerance threshold. 
#' @param genotype_df atav genotypes data frame
#' @param sample_df atav peds file 
#' @param analysis_label_var text string indicating label of analysis. Used to mark saved permutations and load them back in 
#' @param permutation_directory directory for permutations
#' @param output_str save directory for figures
#' @param nperm number of permutations to either create or to load in for display figure
#' @param save_stat determines the type of statistic saved as well as type of analysis performed. Options are "or" where just the OR is saved. "stat" means that whole comparison including or, p value, and test statistic are all saved.
#' @param test_flag determines type of test that you're running. Options are "cmh" "fet" "log"
#' @param decile_label_var i think location on plot of decile labels. 
#' @param max_y_disp maximum height of plot
#' @param y_break_marks location of labels
#' @return saves plots 
#' @examples
#' figure_gene_loeuf_intolerance_threshold(genotype_df,sample_df, paste0("neuro_refractory_loeuf_UrPTV_min_case_",min_case,"_resolution_",resolution_var),permutation_directory, output_str=output_str,nperm=1000,save_stat = "stat", test_flag = "cmh")
#'  figure_gene_intolerance_threshold(genotype_df,sample_df,paste0("neuro_refractory_loeuf_UrPTV_min_case_",min_case,"_resolution_",resolution_var),permutation_directory,output_str=output_str,nperm=100,perform_permutations = FALSE,save_stat = "stat_p", test_flag = "cmh", decile_label_var = 3.5,max_y_disp = 6.1)
figure_gene_loeuf_intolerance_threshold_repo<-function(genotype_df, sample_df, analysis_label_var,permutation_directory, output_str, nperm, numCores, threshold_var,
                                            save_stat = "stat" ,test_flag=NULL, decile_label_var = 5.5, max_y_disp = 10, y_break_marks = c(0,2.5,5)){
  logr::log_print(paste0("Saving permutations into ",permutation_directory))
  if(is.null(test_flag)){
    logr::log_print("Exit: need to define test flag")
    stop()
  } 
  # ****************
  # initialize variables
  # ****************
  # str(genotype_df)
  if(threshold_var=="LOEUF"){
    loeuf_seq <- unique(sort(genotype_df$gnomAD.oe_lof_upper))
    logr::log_print(sprintf("Looking at %i genes with %i unique LOEUF scores", length(unique(genotype_df$Gene.Name)), number_of_bins <- length(loeuf_seq)))
    gene_list <- mclapply(1:number_of_bins,
                          function(x) return(unique((genotype_df %>% filter(gnomAD.oe_lof_upper <= loeuf_seq[x]))$Gene.Name)),mc.cores = numCores )
    names(gene_list)<-rep(sprintf("Less than or equal to %.3f",loeuf_seq),1)
    
    logr::log_print(analysis_label_var <- paste0(analysis_label_var))
  } else if(threshold_var=="mis_z"){
    loeuf_seq <- sort(unique(round(genotype_df$gnomAD.mis_z * 500) / 500), decreasing = TRUE)
    logr::log_print(sprintf("Looking at %i genes with %i unique mis z scores", length(unique(genotype_df$Gene.Name)), number_of_bins <- length(loeuf_seq)))
    gene_list <- mclapply(1:number_of_bins,
                          function(x) return(unique((genotype_df %>% filter(gnomAD.mis_z >= loeuf_seq[x]))$Gene.Name)),mc.cores = numCores )
    names(gene_list)<-rep(sprintf("Greater than or equal to %.3f",loeuf_seq),1)
    
  }  else {
    error("bad thershold_var value")
    return()
  }
  generate_unbiased_threshold_analysis(genotype_df, sample_df, analysis_label_var,permutation_directory, output_str, nperm, gene_list, loeuf_seq,
                                       save_stat, test_flag, decile_label_var, max_y_disp, y_break_marks, numCores, plot_type = threshold_var)
  
}


figure_gene_loeuf_window<-function(genotype_df, sample_df, analysis_label_var,permutation_directory, output_str, nperm = 100000, size_of_bin,
                                                  save_stat = "stat", test_flag=NULL, decile_label_var = 5.5, max_y_disp = 10, y_break_marks = c(0,2.5,5)){
  logr::log_open(paste0(output_str,"/", 
                        paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_loeuf_permutation_analysis_logfile_",analysis_label_var, "_"),"logfile.log"))
  logr::log_print(paste0("Saving permutations into ",permutation_directory))
  if(is.null(test_flag)){
    logr::log_print("Exit: need to define test flag")
    stop()
  } 
  # ****************
  # initialize variables
  # ****************
  gene_bin_list <- return_gene_list_for_loeuf_win(genotype_df, size_of_bin)
  gene_list <- gene_bin_list$gene_list
  x_val <- gene_bin_list$x_val
  
  
  logr::log_print(analysis_label_var <- paste0(analysis_label_var,"_loeuf_win_bin_size_", size_of_bin))
  
  generate_unbiased_threshold_analysis(genotype_df, sample_df, analysis_label_var,permutation_directory, output_str, nperm, gene_list, x_val,
                                       save_stat, test_flag, decile_label_var, max_y_disp, y_break_marks, numCores, plot_type = "LOEUF")
  
}



return_gene_list_for_loeuf_win <- function(genotype_df, size_of_bin, numCores){
  nrow(genotype_unique_sorted_df <- (genotype_df[order(genotype_df$gnomAD.oe_lof_upper),] %>% distinct(Sample.Name, Gene.Name, .keep_all = TRUE) %>% distinct(Sample.Name, .keep_all = TRUE)))
  end_index <- nrow(genotype_unique_sorted_df) - size_of_bin+1
  gene_list <- mclapply(1:end_index,function(x) genotype_unique_sorted_df$Gene.Name[x:(x+size_of_bin-1)], mc.cores = numCores)
  x_val <- sapply(1:end_index,function(x) genotype_unique_sorted_df$gnomAD.oe_lof_upper[round((x+size_of_bin/2))])
  names(gene_list) <- sprintf("Bin centered on LOEUF score %.3f",x_val)
  gene_bin_list <- list(gene_list,x_val)
  names(gene_bin_list) <- c("gene_list","x_val")
  return(gene_bin_list)
}

#' generates analysis
#' @param genotype_df atav genotypes data frame
#' @param sample_df atav peds file 
#' @param analysis_label_var text string indicating label of analysis. Used to mark saved permutations and load them back in 
#' @param permutation_directory directory for permutations
#' @param output_str save directory for figures
#' @param nperm number of permutations to either create or to load in for display figure
#' @param x_val the x axis axis of the plot. corresponds to threshold value of loeuf score, misz, mtr, etc.
#' @param save_stat determines the type of statistic saved as well as type of analysis performed. Options are "or" where just the OR is saved. "stat" means that whole comparison including or, p value, and test statistic are all saved.
#' @param test_flag determines type of test that you're running. Options are "cmh" "fet" "log"
#' @param decile_label_var i think location on plot of decile labels. 
#' @param max_y_disp maximum height of plot
#' @param y_break_marks location of labels#' @return saves plots 
#' @examples 
#' generate_unbiased_threshold_analysis <- function(genotype_df, sample_df, analysis_label_var,permutation_directory, output_str, nperm, gene_list, x_val, save_stat, test_flag, decile_label_var, max_y_disp, y_break_marks)
generate_unbiased_threshold_analysis <- function(genotype_df, sample_df, analysis_label_var,permutation_directory, output_str, nperm, gene_list, x_val,
                                                 save_stat, test_flag, decile_label_var, max_y_disp, y_break_marks, numCores, plot_type = "base"){
  # ****************
  # determine if permutations are needed and then run
  # ****************
  run_perms_if_needed(genotype_df, sample_df, analysis_label_var,permutation_directory, nperm, test_flag, gene_list, numCores)
  if(plot_type == "none"){
    logr::log_print("Returning after finishing permutations. No plots")
    return()
  }
    
  # ****************
  # obtain actual OR vector. I have to this twice because I need the OR and the p vals?
  # ****************
  logr::log_print("Getting non-permuted data")
  actual_matrix_list <- print_or_from_gene_sets(genotype_df, sample_df, gene_list, save_stat, test_flag, nperm,x_val,output_str,analysis_label_var,numCores)
  
  
  actual_or_matrix <- actual_matrix_list$actual_or_matrix
  if(save_stat %in% c("stat", "stat_onesided")){
    matrix_to_gen_p_vals <- actual_matrix_list$actual_or_matrix
  } else if(save_stat %in% c("stat_p")){
    matrix_to_gen_p_vals <- actual_matrix_list$actual_p_matrix
  } else if(save_stat %in% c("perc_diff")){
    matrix_to_gen_p_vals <- actual_matrix_list$actual_perc_diff_matrix
  }
  actual_or_df <- actual_matrix_list$actual_or_df

  # ****************
  # load in permutations and create p values
  # ****************
  
  logr::log_print(sprintf("%s is the name of the saved permutation summary",rds_name <- paste0(permutation_directory,"threshold_permutation_",analysis_label_var,"_nperm_",as.character(nperm),"_",save_stat,".RDS")))
  if(file.exists(rds_name)){
    logr::log_print("Permutation list does exist")
    
    # new_p_vals <- sapply(1:number_of_bins, function(x) return(sum(actual_p_matrix[,x] > permuted_p[,x]))) / nperm
    p_list <- readRDS(file = rds_name)
    logr::log_print("Permutation list loaded")
  } else {
    logr::log_print("Permutation list doesnt exist, running now")
    permuted_stat_matrix <- load_in_gene_list_permutations(nperm,permutation_directory, analysis_label_var, save_stat, length(x_val), numCores)
    p_list <- return_p_val_display_df_list(matrix_to_gen_p_vals, permuted_stat_matrix,nperm,save_stat)
    saveRDS(object = p_list, file = rds_name)
    logr::log_print("Permutation list saved")
  }
  print("debug")
  # min_p <- min(p_list$new_p_vals_vector)
  # min_p_index <- which(p_list$new_p_vals_vector == min_p[1])
  # logr::log_print(sprintf("Smallest permuted numerators is %i which is an unadjusted p value of %.1e. If <= 1, probably need more permutations",
  #                         min_p[1],p_list$permutation_p_vals[min_p_index[1]]))
  # if(length(min_p_index) > 1){
  #   logr::log_print("Need more permutations")
  # }  else {
  #   logr::log_print("Do not need more permutations")
  # }
  # 
  actual_or_df$permutation_p_vals <- p_list$permutation_p_vals
  actual_or_df$new_p_vals_vector <- p_list$new_p_vals_vector
  actual_or_df$perc_diff <- actual_matrix_list$actual_perc_diff_matrix[1,]
  actual_or_df$x_val <- x_val
  actual_or_df_high_quality <- actual_or_df %>% filter((as.numeric(`Ctrl w/ QV`) +  as.numeric(`Case w/ QV`))  <= (1 * nrow(sample_df)))
                                                         # ) #[as.numeric(actual_or_df$`Ctrl w/ QV`) < as.numeric(actual_or_df$`Ctrl w/o QV`),]
  actual_or_df_high_quality$adjusted_p <- p.adjust(actual_or_df_high_quality$permutation_p_vals,method = "fdr")
  actual_or_df_high_quality$adjusted_p_neg_log <- -1*log10(actual_or_df_high_quality$adjusted_p)
  actual_or_df_high_quality$adjusted_p_neg_log[actual_or_df_high_quality$adjusted_p_neg_log==Inf] <- 9
  actual_or_df_high_quality$unadjusted_p_neg_log <- -1*log10(actual_or_df_high_quality$permutation_p_vals)
  actual_or_df_high_quality$unadjusted_p_neg_log[actual_or_df_high_quality$unadjusted_p_neg_log==Inf] <- 9
  save_path_root <- paste0(output_str,plot_type,"_permutation_",analysis_label_var,"_nperm_",as.character(nperm),"_",save_stat)
  
  print("debug")
  min_p <- min(actual_or_df_high_quality$new_p_vals_vector)
  min_p_index <- which(actual_or_df_high_quality$new_p_vals_vector == min_p[1])
  logr::log_print(sprintf("%i minimum P values.",length(min_p_index)))
  logr::log_print(sprintf("The first minimum p-value occurs at index %i which is at x value %0.3f", min_p_index[1], actual_or_df_high_quality$x_val[min_p_index[1]]))
  logr::log_print(sprintf("Smallest permuted numerators is %i which is an unadjusted p value of %.1e. If <= 1, probably need more permutations",
                          min_p[1],actual_or_df_high_quality$permutation_p_vals[min_p_index[1]]))

  logr::log_print(sprintf("The last minimum p-value occurs at index %i which is at x value %0.3f", min_p_index[length(min_p_index)], actual_or_df_high_quality$x_val[min_p_index[length(min_p_index)]]))
  
  if(length(min_p_index) > 1){
    logr::log_print("Need more permutations")
  }  else {
    logr::log_print("Do not need more permutations")
  }
  
  
  
  plot_threshold_output(genotype_df,actual_or_df_high_quality,save_path_root, decile_label_var, max_y_disp, y_break_marks, plot_type = plot_type)
  
  # display_df <- data.frame(x_val = as.numeric(as.character(x_val)), permuted_p_vals = p_list$  , y_val = (-1*log10(new_p_vals)))
  # 
  # min_p <- min(new_p_vals)
  # min_p_loeuf_score <- x_val_var[min_p_index[1]]
  # 
  # display_df$y_val <- gsub(pattern = Inf,replacement =  subst_p_val,x = (-1*log10(new_p_vals)))
  # # fdr
  # display_df$adjusted_p <- p.adjust(new_p_vals,method = p_correction)
  # display_df$adjusted_p_neg_log <- -1*log10(p.adjust(new_p_vals,method = p_correction))
  # display_df$y_val2 <- gsub(pattern = Inf,replacement =  subst_p_val,x = (-1*log10(display_df$adjusted_p)))
  # display_df$actual_or <- actual_or_matrix[1,min_index:max_index]
  # 
  # p_display_list$display_df$p_vals_not_permuted <- actual_matrix_list$actual_p_matrix[1,]
  # debug
  # sum(as.numeric(as.character(p_display_list$display_df$p_vals_not_permuted[752])) > as.numeric(as.character(p_display_list$permuted_or_matrix[,752])))
  # 5/2000
  
}

#'
run_perms_if_needed <- function(genotype_df, sample_df, analysis_label_var,permutation_directory, nperm, test_flag, gene_list, numCores){
  nperm_vector_base <- 1:nperm
  logr::log_print(paste0("Permutation directory is ",permutation_directory))
  temp <- return_files_for_threshold(permutation_directory, analysis_label_var, "stat")
  # extract the permutations which exist of the ones queries
  success <- temp[[2]]
  # determine if any required permutations are not present
  nperm_vector <- nperm_vector_base[!(nperm_vector_base %in% success)]
  
  if(length(nperm_vector)>0){
    logr::log_print("Running remaining permutations starting with these permutations...")
    logr::log_print(head(nperm_vector))
    run_gene_list_permutations(nperm_vector,genotype_df,sample_df,gene_list,analysis_label_var, permutation_directory,"stat",test_flag, numCores)
  }  else {
    logr::log_print("all permutations are complete")
  }
}


#' Determines what permutations have already been done and returns paths to those
#' @param permutation_directory directory for permutations
#' @param analysis_label_var text string indicating label of analysis. Used to mark saved permutations and load them back in 
#' @param save_stat determines the type of statistic saved as well as type of analysis performed. Options are "or" where just the OR is saved. "stat" means that whole comparison including or, p value, and test statistic are all saved.
#' Of note, there is a global variable called mac_pc which is used to determine whether the code is being run on the server or being run locally with mapped folders. This is the way Josh runs the code and I define this variable in the executable script. Not sure if this is the way to go moving forward.
return_files_for_threshold<-function(permutation_directory, analysis_label_var, save_stat){
  logr::log_print("Returning list of file locations of permutations to be loaded in...")
  success <- NULL
  if(save_stat %in% c("stat","stat_p","stat_onesided")){
    files <- list.files(paste0(permutation_directory), pattern = glob2rx(paste0(analysis_label_var, "_perm_", "*","_stat.RDS")))
    success <- as.numeric(gsub(paste0(analysis_label_var, "_perm_|_stat.RDS"), "", files))
  } else {
    logr::log_print(sprintf("save_stat is %s which is not defined", save_stat))# should probably 
    stop()
  }
  
  # if(mac_pc == "server"){  
    full_files <- paste0(permutation_directory,files)
    x <- file.info(full_files)
    print(head(x[order(-x$size),] %>% select(size)))
  # }
  return(list(files,success))
}

#' starts running the permutations. 
#' @param nperm_vector 
#' @param genotype_df 
#' @param sample_df 
#' @param gene_list 
#' @param analysis_label_var 
#' @param permutation_directory 
#' @param save_stat 
#' @param test_flag 
run_gene_list_permutations<-function(nperm_vector,genotype_df,sample_df,gene_list,analysis_label_var, permutation_directory,save_stat, test_flag, numCores){
  
  logr::log_print("Starting Permutations")
  logr::log_print(nperm_vector[nperm_index <- 1:(3*numCores)])
  while(nperm_index[numCores] <= nperm_vector[length(nperm_vector)]){
    mclapply(nperm_vector[nperm_index],
             loeuf_thresholding,
             genotype_df = genotype_df,
             sample_df = sample_df,
             gene_list = gene_list,
             analysis_label = analysis_label_var,
             output_str = permutation_directory,
             permute_var = TRUE,  save_stat=save_stat,
             test_flag = test_flag,
             mc.cores = numCores,mc.preschedule = FALSE)
    
    print(nperm_vector[nperm_index[1]])
    nperm_index <- nperm_index + 3*numCores
  }
  
  # loeuf_thresholding(nperm_vector[1], genotype_df = genotype_df, sample_df = sample_df,   gene_list = gene_list,   analysis_label = analysis_label_var, output_str = permutation_directory,  permute_var = TRUE,  save_stat=save_stat,  test_flag = test_flag)
}

#' Function which thresholds. THIS FUNCTION WILL NOT WORK WITHOUT OTHER FUNCTIONS THAT JOSH HAS CREATED. NEEDS BETTER FACTORING.
#' @param dummy_var dummy variable for permutations
#' @param genotype_df atav genotypes data frame
#' @param sample_df atav peds file 
#' @param gene_list
#' @param output_str
#' @param analysis_label
#' @param permute_var 
#' @param save_stat 
loeuf_thresholding <- function(dummy_var, genotype_df, sample_df, gene_list, output_str = NULL, analysis_label = NULL, permute_var = FALSE, save_stat = "or",
                               cores = 1, test_flag = NULL){
  if(is.null(test_flag)) stop("need to define test_flag")
  if(permute_var == TRUE){
    permuted_dfs <- return_shuffled_atav_cluster_fp(1, atav_temp = genotype_df, sample_temp = sample_df)
    genotype_df <- permuted_dfs$genotype_df %>% select(Variant.ID, Sample.Name, Sample.Phenotype, Gene.Name, Effect, cluster)
    sample_df <- permuted_dfs$sample_df
    if(dummy_var %% 5 == 1) print(dummy_var)
  }
  
  number_of_bins <- length(gene_list)
  loo_cluster_af=rep(NA, number_of_bins)
  
  # using effects that aren't missense or LOF here but would consider tossing
  condition_vector <-rep(c("Unfiltered"),number_of_bins)
  names(condition_vector) <- rep(c( "Unfiltered"),number_of_bins)#
  pub_names <- rep(c("Unfiltered"),number_of_bins)
  
  title_str <-sprintf("%i division",number_of_bins)
  pub_title <- "experiment"
  display_df <- figure_or_forest_fp(genotype_df,sample_df, rep("case",number_of_bins), rep("ctrl",number_of_bins),gene_list,condition_vector,
                                 pub_title,pub_names=pub_names, cores = cores,
                                 loo_cluster_af=loo_cluster_af,mc_meth="none",test_flag = test_flag ) %>% map_df(rev)

  if(permute_var == TRUE){
   if(save_stat %in% c("stat","dec_stat","stat_p","perc_diff")){
      saveRDS(display_df, file = paste0(output_str,analysis_label,"_perm_",as.character(dummy_var),"_",save_stat,".RDS"))
      return((display_df))
    } else {
      logr::log_print("Bad save_stat")
      stop()
    }
  } else{
    return(display_df)
  }
}

return_perc_diff <- function(actual_or_df){
  return(as.numeric(actual_or_df$`Case w/ QV`)/(as.numeric(actual_or_df$`Case w/ QV`) + as.numeric(actual_or_df$`Case w/o QV`)) - 
    as.numeric(actual_or_df$`Ctrl w/ QV`)/(as.numeric(actual_or_df$`Ctrl w/ QV`) + as.numeric(actual_or_df$`Ctrl w/o QV`)))
  
}

#' 
print_or_from_gene_sets <- function(genotype_df, sample_df, gene_list, save_stat, test_flag, nperm,x_val_var,output_str,analysis_label_var, numCores, plot_val = TRUE){
  logr::log_print(sprintf("Starting actual OR with save_stat: %s, and test_flag: %s, and nperm: %i", save_stat, test_flag, nperm))
  actual_or_df <- loeuf_thresholding(1,genotype_df,sample_df, gene_list, permute_var = FALSE,save_stat=save_stat, test_flag = test_flag, cores = numCores) 
  
  actual_or <-as.numeric(as.character(actual_or_df$OR))
  actual_or_matrix <- matrix(rep(actual_or,nperm),nrow=nperm,ncol = length(x_val_var), byrow=TRUE)
  actual_perc_diff <-abs(return_perc_diff(actual_or_df))
  actual_perc_diff_matrix <- matrix(rep(actual_perc_diff,nperm),nrow=nperm,ncol = length(x_val_var), byrow=TRUE)
  actual_p <-as.numeric(as.character(actual_or_df$P_uncorrected))
  actual_p_matrix <- matrix(rep(actual_p,nperm),nrow=nperm,ncol = length(x_val_var), byrow=TRUE)
  
  if(plot_val){
    display_df <- data.frame(x_val = x_val_var, or_vals = as.numeric(as.character(actual_or_df$OR)))
    display_df$or_vals <- gsub(pattern = Inf,replacement =  100,x = display_df$or_vals)
    
    png(paste(output_str,"actual_data_",analysis_label_var,"_nperm_",as.character(nperm),"_",save_stat, ".png",sep = ""), width = 1200, height = 1200)
    temp <- ggplot(display_df, aes(as.numeric(as.character(x_val)),as.numeric(as.character(or_vals)))) +
      geom_line() + 
      theme_cowplot(20)
    print(temp)
    dev.off()
  }
  
  
  return_list <- list(actual_or_df, actual_or_matrix,actual_p_matrix,actual_perc_diff_matrix)
  names(return_list) <- c("actual_or_df","actual_or_matrix","actual_p_matrix","actual_perc_diff_matrix")
  
  return(return_list)
  
}

#' 
return_p_val_display_df_list <- function(actual_stat_matrix, permuted_stat_matrix,nperm,save_stat){
  # permuted_or_matrix_temp <- gsub(0,Inf,permuted_or_matrix)
  # inf_index <- which(actual_or_matrix[1,] == Inf)
  # ifelse(length(inf_index) > 0, max_index <- inf_index[1] - 1,max_index <- length(actual_or_matrix[1,]))
  # if(!is.null(actual_or_df)){
  #   onehundred_index <- min(which(actual_or_df$`Case w/o QV` == 0),which(actual_or_df$`Ctrl w/o QV` == 0))
  #   max_index <- min(onehundred_index, max_index)
  #   logr::log_print(sprintf("max index accounting for 100 percent saturated is %i",max_index))
  # }
  # logr::log_print(sprintf("Length of vector is %i, max index accounting for infinity is %i", length(actual_or_matrix[1,] ),max_index))
  # if(!is.null(actual_or_df)){
  #   fifty_percent <- min(which(actual_or_df$`Ctrl % w/ QV` >50))
  #   max_index <- min(fifty_percent, max_index)
  #   logr::log_print(sprintf("max index accounting for controls are 50% saturated is %i",max_index))
  # }
  
  min_index <- 1
  if(save_stat %in% c("or","stat")){
    # this is a two sided test, either with odds ratio or i think the test statistic?
    new_p_vals_vector <- sapply(1:length(actual_stat_matrix[1,]), function(x) return(sum(ifelse(actual_stat_matrix[,x] > 1 ,actual_stat_matrix[,x],1/actual_stat_matrix[,x] ) <= ifelse(permuted_stat_matrix[,x] > 1 ,permuted_stat_matrix[,x],1/permuted_stat_matrix[,x] ))))
    # new_p_vals_vector <- sapply(1:max_index, function(x) return(sum(actual_or_matrix[,x] <= permuted_or_matrix[,x])))
  } else if (save_stat %in% c("stat_onesided","perc_diff")){
    # this is a one sided test, i think it gets used with the odds ratio?
    new_p_vals_vector <- sapply(1:length(actual_stat_matrix[1,]), function(x) return(sum(actual_stat_matrix[,x] <= permuted_stat_matrix[,x]))) 
  } else if (save_stat == "stat_p"){
    # This is the most common method of running this code. Here, we calculate a one sided permutation test where the permuted values are p values. I am calculating how often the actual p value is greater than the permuted p values. That is my permuted P value.
    new_p_vals_vector <- sapply(1:length(actual_stat_matrix[1,]), function(x) return(sum(actual_stat_matrix[,x] >= permuted_stat_matrix[,x])))
  } 
  
  permutation_p_vals <-  new_p_vals_vector/ nperm

  return_list <- list(permutation_p_vals,actual_stat_matrix, permuted_stat_matrix, save_stat,nperm, new_p_vals_vector)
  names(return_list) <- c("permutation_p_vals","actual_stat_matrix", "permuted_stat_matrix", "save_stat","nperm","new_p_vals_vector")
  return(return_list)
  
}

#' This function loads in all of the permutations. 
load_in_gene_list_permutations<-function(nperm,permutation_directory, analysis_label_var, save_stat, number_of_bins, numCores){
  return_permuted_stat <- function(files){
    temp_stat <- mclapply(files, readRDS, mc.cores = numCores)
    if(save_stat %in% c("stat","stat_onesided")){
      permuted_stat_pre_inf <- as.numeric(sapply(1:length(temp_stat),function(x) temp_stat[[x]]$OR),simplify = TRUE)#, byrow=TRUE)
    } else if(save_stat == "stat_p"){
      permuted_stat_pre_inf <- as.numeric(sapply(1:length(temp_stat),function(x) temp_stat[[x]]$P_uncorrected),simplify = TRUE)#, byrow=TRUE)
      # permuted_or[permuted_or==Inf] <- 10000000
      # josh_debug <- as.numeric(sapply(1:length(temp_stat),function(x) temp_stat[[x]]$P_uncorrected[500]),simplify = TRUE)#, byrow=TRUE)
    } else if(save_stat == "perc_diff"){
      permuted_stat_pre_inf <- sapply(1:length(temp_stat),function(x) abs(return_perc_diff(temp_stat[[x]])),simplify = TRUE)#, byrow=TRUE)
    }
    permuted_stat <- permuted_stat_pre_inf
    permuted_stat[permuted_stat==Inf] <- 10000000#corrects if OR is infinity
    return(permuted_stat_matrix <- matrix(permuted_stat,nrow=length(files),ncol = number_of_bins, byrow=TRUE))
  }
  
  
  logr::log_print("Loading in saved Permutations in load_in_gene_list_permutations...")
  nperm_vector <- 1:nperm
  files_all <- paste0(permutation_directory, analysis_label_var, "_perm_", 1:nperm,"_stat.RDS")
  # do in 10K aliquots
  denominator_val <- 10000
  num_10k_aliquots <- as.integer(nperm/denominator_val)
  remainder_of_perms <- nperm %% denominator_val
  if(remainder_of_perms != 0){
    permuted_stat_matrix <- return_permuted_stat(files_all)
  } else {
    start_indices <- seq(from = 1, to = nperm, by = denominator_val)
    stop_indices <- seq(from = denominator_val, to = nperm, by = denominator_val)
    if(length(start_indices) != length(stop_indices)){
      stop()
    }
    list_of_permuted_stat_pre_inf <- lapply(1:length(start_indices), function(x) return_permuted_stat(files_all[start_indices[x]:stop_indices[x]]))
    permuted_stat_matrix <- do.call(rbind, list_of_permuted_stat_pre_inf)
  }
  return(permuted_stat_matrix)
}

# Functions for displaying--------
#' Display the threshold data
#' @param genotype_df atav genotypes data frame
#' @param p_display_list 
#' @param x_val 
#' @param actual_or_matrix 
#' @param plot_type determines the type of plot so can change text labels, etc 
plot_threshold_output<-function(genotype_df, display_df,save_path_root, decile_label_var, max_y_disp, y_break_marks, plot_type = "base"){
  display_df$case_ctrl <- ifelse(as.numeric(display_df$OR > 1), "Case","Control")
  
  min_p <- min(display_df$new_p_vals_vector)
  min_p_index <- which(display_df$new_p_vals_vector == min_p[1])
  # min_p <- min(display_df$adjusted_p)
  # min_p_index <- which(display_df$adjusted_p == min_p[1])
  
  x <- 1
  # CI_top <- sapply(1:length(new_p_vals), function(x) sort(p_display_list$permuted_or_matrix[x,])[floor(.05 * length(new_p_vals))])
  # CI_top <- sapply(1:nrow(display_df), function(x) sort(p_display_list$permuted_or_matrix[,x])[floor(.05 * length(new_p_vals))]) #after loading in a permutation set, the dimension s of the permuted_or_matrix seemed to be switched. I swapped the x to the other dimension here. Not sure if this will cause a bug later
  # display_df$CI_top <- CI_top
  
  
  box_colors <- rep(c("gray76","gray33"),5)
  alpha_vector <- rep(0.5,10)
  
  theme_set(theme_cowplot(font_family = "ArialMT"))
  # ***************
  # This plot points to the most sigificant threshold
  # ***************
  adjusted_p_plot_fig <- ggplot(display_df, aes(as.numeric(as.character(x_val)),as.numeric(as.character(adjusted_p_neg_log)))) +
    # geom_line(size = .75) + 
    geom_point(size = .25) + 
    geom_hline(yintercept = -1 * log10(0.05), color = "red", size = .75) + 
    # geom_label(aes(x = 1, y = decile_label_var+2), label = cat("'Test '", parse(text=expression_string))) + #,pretty10exp())) display_df$adjusted_p[min_p_index[1]]
    annotate(geom = "point",x = display_df$x_val[min_p_index[1]], y = as.numeric(display_df$adjusted_p_neg_log[min_p_index[1]]), color = "red", size = 3) +
    ylab(bquote('-Log'[10]~'(P-Value)')) +
    # ylab("-Log-10(p-value)") +
    theme_cowplot(font_size = 12) +
    # gghighlight::gghighlight(x_val == min_p_loeuf_score, ) +
    scale_y_continuous(limits = c(0,max_y_disp),expand = c(0,0), breaks = y_break_marks) +
    scale_x_continuous(expand = c(0,0)) +
    theme(plot.margin = unit(rep(c(0.25),4), "in"))
  
  if(plot_type == "LOEUF"){
    x_axis_label <- "LOEUF"
    decile_max <- sapply(0:9,
                         function(x) return(max((genotype_df %>% filter(gnomAD.oe_lof_upper_bin == x))$gnomAD.oe_lof_upper, na.rm = TRUE)))
    decile_min <- sapply(0:9,
                         function(x) return(min((genotype_df %>% filter(gnomAD.oe_lof_upper_bin == x))$gnomAD.oe_lof_upper, na.rm = TRUE)))
    decile_label_x <- sapply(1:10,function(x) mean(c(decile_max[x],decile_min[x])))
    decile_index <- !(decile_max %in% c(Inf,-Inf)) & !(decile_min %in% c(Inf,-Inf))
    adjusted_p_plot_fig <- adjusted_p_plot_fig + 
      annotate("text",x=decile_label_x,y=decile_label_var, label = 1:10, size = 5) +
      annotate("rect",xmin = decile_min[decile_index],xmax = decile_max[decile_index],ymin = rep(decile_label_var - 0.6,10)[decile_index],
               ymax = rep(decile_label_var - 0.5,10)[decile_index], fill =  box_colors[decile_index], alpha = alpha_vector[decile_index] ) +
      geom_vline(xintercept = decile_max, color = "black", size = 0.5, linetype = 2, alpha = 0.5)
    
    
    
  } else if (plot_type == "mis_z"){
    x_axis_label <- "Missense Z"
    adjusted_p_plot_fig <- adjusted_p_plot_fig + scale_x_reverse()
  } else {
    x_axis_label <- "NO X AXIS DEFINED, PLOT TYPE BASE"
  }
  adjusted_p_plot_fig <- adjusted_p_plot_fig + 
    xlab(x_axis_label) +
    geom_label(aes(x = 1, y = decile_label_var+2, label = sprintf("Most significant OR is %s at %s = %.3f\nwith an adjusted p-value of %.1e",display_df$OR[min_p_index[1]],x_axis_label,x_val[min_p_index[1]],display_df$adjusted_p[min_p_index[1]])))
  
  ggsave(plot = adjusted_p_plot_fig, filename = paste(save_path_root, "_most_sig_threshold_adjusted_p.pdf",sep = ""), width = 4.488, height = 3.5, units = "in", dpi = 300)
  
  # ***************
  # Same Plot but with - log of raw p values
  # ***************unadjusted_p_neg_log
  unadjusted_p_plot_fig <- ggplot(display_df, aes(as.numeric(as.character(x_val)),as.numeric(as.character(unadjusted_p_neg_log)))) +
    # geom_line(size = .75) + 
    geom_point(size = .25) + 
    geom_hline(yintercept = -1 * log10(0.05/nrow(display_df)), color = "red", size = .75) + 
    # geom_label(aes(x = 1, y = decile_label_var+2), label = cat("'Test '", parse(text=expression_string))) + #,pretty10exp())) display_df$adjusted_p[min_p_index[1]]
    annotate(geom = "point",x = display_df$x_val[min_p_index[1]], y = as.numeric(as.character(display_df$unadjusted_p_neg_log[min_p_index[1]])), color = "red", size = 3) +
    ylab(bquote('-Log'[10]~'(P-Value)')) +
    # ylab("-Log-10(p-value)") +
    theme_cowplot(font_size = 12) +
    # gghighlight::gghighlight(x_val == min_p_loeuf_score, ) +
    scale_y_continuous(limits = c(0,max_y_disp),expand = c(0,0), breaks = y_break_marks) +
    scale_x_continuous(expand = c(0,0)) +
    theme(plot.margin = unit(rep(c(0.25),4), "in"))
  
  if(plot_type == "LOEUF"){
    x_axis_label <- "LOEUF"
    decile_max <- sapply(0:9,
                         function(x) return(max((genotype_df %>% filter(gnomAD.oe_lof_upper_bin == x))$gnomAD.oe_lof_upper, na.rm = TRUE)))
    decile_min <- sapply(0:9,
                         function(x) return(min((genotype_df %>% filter(gnomAD.oe_lof_upper_bin == x))$gnomAD.oe_lof_upper, na.rm = TRUE)))
    decile_label_x <- sapply(1:10,function(x) mean(c(decile_max[x],decile_min[x])))
    decile_index <- !(decile_max %in% c(Inf,-Inf)) & !(decile_min %in% c(Inf,-Inf))
    unadjusted_p_plot_fig <- unadjusted_p_plot_fig + 
      annotate("text",x=decile_label_x,y=decile_label_var, label = 1:10, size = 5) +
      annotate("rect",xmin = decile_min[decile_index],xmax = decile_max[decile_index],ymin = rep(decile_label_var - 0.6,10)[decile_index],
               ymax = rep(decile_label_var - 0.5,10)[decile_index], fill =  box_colors[decile_index], alpha = alpha_vector[decile_index] ) +
      geom_vline(xintercept = decile_max, color = "black", size = 0.5, linetype = 2, alpha = 0.5)
    
    
    
  } else if (plot_type == "mis_z"){
    x_axis_label <- "Missense Z"
    unadjusted_p_plot_fig <- unadjusted_p_plot_fig + scale_x_reverse()
  } else {
    x_axis_label <- "NO X AXIS DEFINED, PLOT TYPE BASE"
  }
  unadjusted_p_plot_fig <- unadjusted_p_plot_fig + 
    xlab(x_axis_label) +
    geom_label(aes(x = 1, y = decile_label_var+2, label = sprintf("Most significant OR is %s at %s = %.3f\nwith an unadjusted p-value of %.1e",display_df$OR[min_p_index[1]],x_axis_label,x_val[min_p_index[1]],display_df$permutation_p_vals[min_p_index[1]])))
  
  ggsave(plot = unadjusted_p_plot_fig, filename = paste(save_path_root, "_most_sig_threshold_unadjusted_p.pdf",sep = ""), width = 4.488, height = 3.5, units = "in", dpi = 300)
  
  return()
  
  # ***************
  # This plot points to the most intolerant threshold that is significant
  # ***************
  display_df$index <- 1:nrow(display_df)
  all_sig_p_vals <- display_df$x_val[display_df$adjusted_p < 0.05]
  min_sig_p_index <- display_df$index[display_df$x_val == all_sig_p_vals[1]]
  temp <- ggplot(display_df, aes(as.numeric(as.character(x_val)),as.numeric(as.character(y_val2)))) +
    # geom_line(size = .75) + 
    geom_point(size = .25) + 
    geom_hline(yintercept = -1 * log10(0.05), color = "red", size = .75) + 
    # geom_label(aes(x = 1, y = decile_label_var+2), label = cat("'Test '", parse(text=expression_string))) + #,pretty10exp())) display_df$adjusted_p[min_p_index[1]]
    annotate(geom = "point",x = display_df$x_val[min_sig_p_index], y = as.numeric(display_df$y_val2[min_sig_p_index]), color = "red", size = 3) +
    ylab(bquote('-Log'[10]~'(P-Value)')) +
    # ylab("-Log-10(p-value)") +
    theme_cowplot(font_size = 12) +
    # gghighlight::gghighlight(x_val == min_p_loeuf_score, ) +
    scale_y_continuous(limits = c(0,max_y_disp),expand = c(0,0), breaks = y_break_marks) +
    scale_x_continuous(expand = c(0,0)) +
    theme(plot.margin = unit(rep(c(0.25),4), "in"))
  
  if(plot_type == "LOEUF"){
    x_axis_label <- "LOEUF"
    decile_max <- sapply(0:9,
                         function(x) return(max((genotype_df %>% filter(gnomAD.oe_lof_upper_bin == x))$gnomAD.oe_lof_upper, na.rm = TRUE)))
    decile_min <- sapply(0:9,
                         function(x) return(min((genotype_df %>% filter(gnomAD.oe_lof_upper_bin == x))$gnomAD.oe_lof_upper, na.rm = TRUE)))
    decile_label_x <- sapply(1:10,function(x) mean(c(decile_max[x],decile_min[x])))
    decile_index <- !(decile_max %in% c(Inf,-Inf)) & !(decile_min %in% c(Inf,-Inf))
    temp <- temp + 
      annotate("text",x=decile_label_x,y=decile_label_var, label = 1:10, size = 5) +
      annotate("rect",xmin = decile_min[decile_index],xmax = decile_max[decile_index],ymin = rep(decile_label_var - 0.6,10)[decile_index],
               ymax = rep(decile_label_var - 0.5,10)[decile_index], fill =  box_colors[decile_index], alpha = alpha_vector[decile_index] ) +
      geom_vline(xintercept = decile_max, color = "black", size = 0.5, linetype = 2, alpha = 0.5)
    
    
    
  } else if (plot_type == "misz"){
    x_axis_label <- "Missense Z"
  } else {
    x_axis_label <- "NO X AXIS DEFINED, PLOT TYPE BASE"
  }
  temp <- temp + 
    xlab(x_axis_label) +
    geom_label(aes(x = 1, y = decile_label_var+2, label = sprintf("The most intolerant threshold\nwith a significant OR is %.1f at %s = %.3f\nwith an adjusted p-value of %.1e",actual_or_matrix[1,min_sig_p_index],x_axis_label,x_val[min_sig_p_index],display_df$adjusted_p[min_sig_p_index])))
  
  ggsave(plot = temp, filename = paste(save_path_root, "_most_intolerant_threshold.pdf",sep = ""), width = 4.488, height = 3.5, units = "in", dpi = 300)  
  
  
  # plot actual P values
  display_df <- p_display_list$display_df
  # ***************
  # This plot points uses the actual p values but applies 
  # ***************
  display_df$adjusted_p <- p.adjust(p_display_list$p_vals_not_permuted,method = "fdr")
  display_df$adjusted_p_neg_log <- -1*log10(p.adjust(display_df$adjusted_p,method = "fdr"))
  display_df$y_val2 <- display_df$adjusted_p_neg_log
  # display_df$y_val2 <- gsub(pattern = Inf,replacement =  subst_p_val,x = (-1*log10(display_df$adjusted_p)))
  display_df$index <- 1:nrow(display_df)
  all_sig_p_vals <- display_df$x_val[display_df$adjusted_p < 0.05]
  min_sig_p_index <- display_df$index[display_df$x_val == all_sig_p_vals[1]]
  temp <- ggplot(display_df, aes(as.numeric(as.character(x_val)),as.numeric(as.character(y_val2)))) +
    # geom_line(size = .75) + 
    geom_point(size = .25) + 
    geom_hline(yintercept = -1 * log10(0.05), color = "red", size = .75) + 
    # geom_label(aes(x = 1, y = decile_label_var+2), label = cat("'Test '", parse(text=expression_string))) + #,pretty10exp())) display_df$adjusted_p[min_p_index[1]]
    annotate(geom = "point",x = display_df$x_val[min_sig_p_index], y = as.numeric(display_df$y_val2[min_sig_p_index]), color = "red", size = 3) +
    ylab(bquote('-Log'[10]~'(P-Value)')) +
    # ylab("-Log-10(p-value)") +
    theme_cowplot(font_size = 12) +
    # gghighlight::gghighlight(x_val == min_p_loeuf_score, ) +
    scale_y_continuous(limits = c(0,max_y_disp),expand = c(0,0), breaks = y_break_marks) +
    scale_x_continuous(expand = c(0,0)) +
    theme(plot.margin = unit(rep(c(0.25),4), "in"))
  
  if(plot_type == "LOEUF"){
    x_axis_label <- "LOEUF"
    decile_max <- sapply(0:9,
                         function(x) return(max((genotype_df %>% filter(gnomAD.oe_lof_upper_bin == x))$gnomAD.oe_lof_upper, na.rm = TRUE)))
    decile_min <- sapply(0:9,
                         function(x) return(min((genotype_df %>% filter(gnomAD.oe_lof_upper_bin == x))$gnomAD.oe_lof_upper, na.rm = TRUE)))
    decile_label_x <- sapply(1:10,function(x) mean(c(decile_max[x],decile_min[x])))
    decile_index <- !(decile_max %in% c(Inf,-Inf)) & !(decile_min %in% c(Inf,-Inf))
    temp <- temp + 
      annotate("text",x=decile_label_x,y=decile_label_var, label = 1:10, size = 5) +
      annotate("rect",xmin = decile_min[decile_index],xmax = decile_max[decile_index],ymin = rep(decile_label_var - 0.6,10)[decile_index],
               ymax = rep(decile_label_var - 0.5,10)[decile_index], fill =  box_colors[decile_index], alpha = alpha_vector[decile_index] ) +
      geom_vline(xintercept = decile_max, color = "black", size = 0.5, linetype = 2, alpha = 0.5)
    
    
    
  } else if (plot_type == "misz"){
    x_axis_label <- "Missense Z"
  } else {
    x_axis_label <- "NO X AXIS DEFINED, PLOT TYPE BASE"
  }
  temp <- temp + 
    xlab(x_axis_label) +
    geom_label(aes(x = 1, y = decile_label_var+2, label = sprintf("The most intolerant threshold\nwith a significant OR is %.1f at %s = %.3f\nwith an adjusted p-value of %.1e",actual_or_matrix[1,min_sig_p_index],x_axis_label,x_val[min_sig_p_index],display_df$adjusted_p[min_sig_p_index])))
  
  ggsave(plot = temp, filename = paste(save_path_root, "_most_intolerant_threshold_non_permuted_p.pdf",sep = ""), width = 4.488, height = 3.5, units = "in", dpi = 300)  
  
  return()
  # past here are different ways to visualize the data but have not been debugged since code was re-factored.
  
  pdf(paste(save_path_root, "_2.pdf",sep = ""), width = 4.488, height = 3.5, family = font_variable, pointsize = 5)
  temp1 <- gene_list_p_figure(display_df,"LOEUF Score", actual_or_matrix,min_p_index)
  temp2 <- temp1 + 
    # geom_vline(xintercept = decile_max, color = "black", size = 0.5, linetype = 2, alpha = 0.5) + 
    annotate("segment",x = decile_max,xend = decile_max,y = rep(0,10),
             yend = rep(1.3,10)[decile_index], alpha = rep(0.25,length(decile_max)), linetype = 2 ) +
    annotate("text",x=decile_label_x,y=1.5, label = 1:10, size = 5) +
    annotate("rect",xmin = decile_min[decile_index],xmax = decile_max[decile_index],ymin = rep(1.3,10)[decile_index],
             ymax = rep(1.35,10)[decile_index], fill =  box_colors[decile_index], alpha = alpha_vector[decile_index] )
  print(temp2)
  dev.off()
  
  # plot 3
  # print("debug")
  pdf(paste(save_path_root, "_3.pdf",sep = ""), width = 4.488, height = 3.5, family = font_variable, pointsize = 5)
  temp1 <- gene_list_p_figure(display_df,"LOEUF Score", actual_or_matrix,min_p_index,p_var = "unadjusted" )
  temp2 <- temp1 + 
    # geom_vline(xintercept = decile_max, color = "black", size = 0.5, linetype = 2, alpha = 0.5) + 
    annotate("segment",x = decile_max,xend = decile_max,y = rep(0,10),
             yend = rep(1.3,10)[decile_index], alpha = rep(0.25,length(decile_max)), linetype = 2 ) +
    annotate("text",x=decile_label_x,y=1.5, label = 1:10, size = 5) +
    annotate("rect",xmin = decile_min[decile_index],xmax = decile_max[decile_index],ymin = rep(1.3,10)[decile_index],
             ymax = rep(1.35,10)[decile_index], fill =  box_colors[decile_index], alpha = alpha_vector[decile_index] )
  print(temp2)
  dev.off()
  # plot 4
  temp <- ggplot(display_df, aes(as.numeric(as.character(x_val)),as.numeric(as.character(y_val)))) +
    # geom_line(size = .75) + 
    geom_point(size = .25) + 
    geom_point(aes(y = (-1 * log10(CI_top)) ), color = "green", size = 0.1) + 
    annotate("text",x=decile_label_x,y=decile_label_var, label = 1:10, size = 5) +
    annotate("rect",xmin = decile_min[decile_index],xmax = decile_max[decile_index],ymin = rep(decile_label_var - 0.6,10)[decile_index],
             ymax = rep(decile_label_var - 0.5,10)[decile_index], fill =  box_colors[decile_index], alpha = alpha_vector[decile_index] ) +
    geom_hline(yintercept = -1 * log10(0.05/nrow(display_df)), color = "red", size = .75) + 
    geom_vline(xintercept = decile_max, color = "black", size = 0.5, linetype = 2, alpha = 0.5) + 
    geom_label(aes(x = 1, y = decile_label_var+2, 
                   label = sprintf("Most significant OR is %.1f at LOUEF = %.3f\nwith an unadjusted p-value of %.1e",actual_or_matrix[1,min_p_index[1]],loeuf_seq[min_p_index[1]],display_df$p_vals[min_p_index[1]]))) +
    # geom_label(aes(x = 1, y = decile_label_var+2), label = cat("'Test '", parse(text=expression_string))) + #,pretty10exp())) display_df$adjusted_p[min_p_index[1]]
    annotate(geom = "point",x = display_df$x_val[min_p_index[1]], y = as.numeric(display_df$y_val[min_p_index[1]]), color = "red", size = 3) +
    xlab("LOEUF Score") +
    ylab(bquote('-Log'[10]~'(P-Value)')) +
    # ylab("-Log-10(p-value)") +
    theme_cowplot(font_size = 12) +
    # gghighlight::gghighlight(x_val == min_p_loeuf_score, ) +
    scale_y_continuous(limits = c(0,max_y_disp),expand = c(0,0), breaks = y_break_marks) +
    scale_x_continuous(expand = c(0,0)) +
    theme(plot.margin = unit(rep(c(0.25),4), "in"))
  ggsave(plot = temp, filename = paste(save_path_root, "_4.pdf",sep = ""), width = 4.488, height = 3.5, units = "in", dpi = 300)
  
}

#Permutations run for a specific gene set-----------

#' Measure optimal p value for loeuf intolerance threshold by gene bin. 
#' @param genotype_df atav genotypes data frame
#' @param sample_df atav peds file 
#' @param analysis_label_var text string indicating label of analysis. Used to mark saved permutations and load them back in 
#' @param permutation_directory directory for permutations
#' @param output_str save directory for figures
#' @param nperm number of permutations to either create or to load in for display figure
#' @param save_stat determines the type of statistic saved as well as type of analysis performed. Options are "or" where just the OR is saved. "stat" means that whole comparison including or, p value, and test statistic are all saved. NOT YET IMPLEMENTED
#' @param test_flag determines type of test that you're running. Options are "cmh" "fet" "log". NOT YET IMPLEMENTED.
#' @param decile_label_var i think location on plot of decile labels. 
#' @param max_y_disp maximum height of plot
#' @param y_break_marks location of labels
#' @return saves plots 
#' @examples
#' figure_gene_loeuf_intolerance_threshold(genotype_df,sample_df, paste0("neuro_refractory_loeuf_UrPTV_min_case_",min_case,"_resolution_",resolution_var),permutation_directory, output_str=output_str,nperm=1000,save_stat = "stat", test_flag = "cmh")
#'  figure_gene_intolerance_threshold(genotype_df,sample_df,paste0("neuro_refractory_loeuf_UrPTV_min_case_",min_case,"_resolution_",resolution_var),permutation_directory,output_str=output_str,nperm=100,perform_permutations = FALSE,save_stat = "stat_p", test_flag = "cmh", decile_label_var = 3.5,max_y_disp = 6.1)
figure_gene_loeuf_intolerance_threshold_by_gene_bin<-function(genotype_df, sample_df, analysis_label_var,permutation_root, output_str, plot_type, nperm = 100000,
                                                  save_stat = "stat", test_flag=NULL, decile_label_var = 5.5, max_y_disp = 10, y_break_marks = c(0,2.5,5),
                                                  size_of_bin = NA, plot_output = TRUE){
  logr::log_open(paste0(output_str,"/", 
                        paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_loeuf_permutation_analysis_logfile_",analysis_label_var, "_"),"logfile.log"))

  
  if(is.null(test_flag)){
    logr::log_print("Exit: need to define test flag")
    stop()
  } 
  # ****************
  # Initialize gene bins
  # ****************
  if(plot_type == "LOEUF_thresh"){
    x_val <- unique(sort(genotype_df$gnomAD.oe_lof_upper))
    logr::log_print(sprintf("Looking at %i genes with %i unique LOEUF scores", length(unique(genotype_df$Gene.Name)), number_of_bins <- length(x_val)))
    gene_list <- mclapply(1:number_of_bins,
                          function(x) return(unique((genotype_df %>% filter(gnomAD.oe_lof_upper <= loeuf_seq[x]))$Gene.Name)),mc.cores = numCores )
    names(gene_list)<-rep(sprintf("Less than or equal to %.3f",x_val),1)
    analysis_label_var <- paste0(analysis_label_var,"_",plot_type)
  } else if (plot_type == "LOEUF_win"){
    nrow(genotype_unique_sorted_df <- (genotype_df[order(genotype_df$gnomAD.oe_lof_upper),] %>% distinct(Sample.Name, Gene.Name, .keep_all = TRUE) %>% distinct(Sample.Name, .keep_all = TRUE)))
    end_index <- nrow(genotype_unique_sorted_df) - size_of_bin+1
    gene_list <- mclapply(1:end_index,function(x) genotype_unique_sorted_df$Gene.Name[x:(x+size_of_bin-1)], mc.cores = numCores)
    x_val <- sapply(1:end_index,function(x) genotype_unique_sorted_df$gnomAD.oe_lof_upper[round((x+size_of_bin/2))])
    names(gene_list) <- sprintf("Bin centered on LOEUF score %.3f",x_val)
    analysis_label_var <- paste0(analysis_label_var,"_",plot_type,"_", size_of_bin)
  }
  logr::log_print(paste0("Updated analysis var is ",analysis_label_var <- paste0(analysis_label_var,"_by_gene_bin")))
  permutation_directory <-paste0(permutation_root,analysis_label_var,"/")
  dir.create(permutation_directory,showWarnings = TRUE)
  # ****************
  # check to see if p values are saved. If not, run permutations
  # ****************
  rdsName <- paste0(permutation_directory, analysis_label_var,"_saved_p_vals.RDS")
  if(file.exists(rdsName)){
    logr::log_print("Permutation list does exist")
    if(plot_output)
      permuted_p_val_list <- readRDS(rdsName)
  } else {
    logr::log_print("Permutation list does not exist. Running permutations.")
    run_gene_list_permutations_by_bin(permutation_directory,analysis_label_var,nperm, gene_list, genotype_df,sample_df)
    # ****************
    # return permuted p vals
    # ****************
    if(plot_output){
      permuted_p_val_list <- return_permutation_p_vals_from_permutation_by_bin(permutation_directory,analysis_label_var,nperm, gene_list, genotype_df,sample_df)
      saveRDS(object = permuted_p_val_list, file = rdsName)
    }
  }
  if(plot_output){
  
    # ****************
    # plot result
    # ****************
    display_df <- data.frame(x_val = x_val,
                             OR = permuted_p_val_list$estimate_true_vector,
                             adjusted_p = p.adjust(permuted_p_val_list$permutation_p_vals, method = "fdr"))
    display_df$adjusted_p_neg_log <- -1*log10(display_df$adjusted_p)
    display_df$adjusted_p_neg_log[display_df$adjusted_p_neg_log==Inf] <- 10
    display_df$case_ctrl <- ifelse(as.numeric(display_df$OR > 1), "Case","Control")
    
    plot_threshold_output(genotype_df, display_df,paste0(output_str,"test"), 5.5, 10,  c(0,2.5,5), plot_type = "LOEUF")
  }
  return(permuted_p_val_list)
}


#' Performs permutations for a specific gene set
#' @param genotype_df atav genotypes data frame
#' @param sample_df atav peds file 
#' @param analysis_label_var text string indicating label of analysis. Used to mark saved permutations and load them back in 
#' @param permutation_directory directory for permutations
#' @param save_stat determines the type of statistic saved as well as type of analysis performed. Options are "or" where just the OR is saved. "stat" means that whole comparison including or, p value, and test statistic are all saved.
#' @param test_flag determines type of test that you're running. Options are "cmh" "fet" "log". THIS OPTION HASN'T BEEN IMPLEMENTED YET
gene_bin_theshold <- function(dummy_var, genotype_df_orig,sample_df_orig,gene_list,bin_index,output_str,permute_var,analysis_label){
  if(permute_var == TRUE){
    # print("in permutation")
    permuted_dfs <- return_shuffled_atav_cluster_fp(1, atav_temp = genotype_df_orig %>% filter(Gene.Name %in% gene_list[[bin_index]]) %>% select(Sample.Name, Sample.Phenotype, Gene.Name, Effect, epi_phenotype, cluster, Variant.ID), sample_temp = sample_df_orig)
    genotype_df <- permuted_dfs$genotype_df
    sample_df <- permuted_dfs$sample_df
  } else{
    # print("no permutation")
    genotype_df <- genotype_df_orig %>% filter(Gene.Name %in% gene_list[[bin_index]]) %>% select(Sample.Name, Sample.Phenotype, Gene.Name, Effect, epi_phenotype, cluster, Variant.ID)
    sample_df <- sample_df_orig
  }
  # print(sample_df)
  stat_output <- atav_fisher(genotype_df, sample_df, c("case","ctrl"),cmh_test=TRUE)
  if(permute_var == TRUE){
    saveRDS(stat_output, file = paste0(output_str,analysis_label,"_bin_", bin_index,"_perm_",as.character(dummy_var),"_stat.RDS"))
    return(stat_output)
  } else{
    return(stat_output)
  }
}

#' Calls function to run permutations for an index of a gene set
run_gene_list_permutations_by_bin <- function(permutation_directory,analysis_label,nperm,gene_list, genotype_df, sample_df, numCores){
  for(bin_index in 1:length(gene_list)){
    print(paste0("Bin Index ", bin_index, " out of ", length(gene_list)))
    success <- NULL
    files <- list.files(paste0(permutation_directory), pattern = glob2rx(paste0(analysis_label,"_bin_", bin_index, "_perm_", "*","_stat.RDS")))
    success <- as.numeric(gsub(paste0(analysis_label, "_bin_", bin_index,"_perm_|_stat.RDS"), "", files))
    print(paste0("Need to do ",length(nperm_vector <- (1:nperm)[1:nperm %!in% success])," out of ", nperm))
    # files <- paste0(permutation_directory,analysis_label,"_bin_", bin_index, "_perm_", 1:nperm,"_stat.RDS")
    # file_exist_var <- sapply(files,file.exists)
    # success <- (1:nperm)[!file_exist_var]
    # temp_stat <- mclapply(files, readRDS, mc.cores = numCores)
    # print("Extracting P values...")
    # p_vals_of_random <- sapply(1:nperm, function(x) temp_stat[[x]][[1]]$p.value)
    # not_permuted <- gene_bin_theshold(1, genotype_df, sample_df,gene_list,bin_index,permutation_directory,FALSE,analysis_label)
    # permutation_p_numerator <- sum(not_permuted[[1]]$p.value >= p_vals_of_random)
    # if(permutation_p_numerator > 2)

    
    if(length(nperm_vector) > 0){
      print("finishing permutations")
      mclapply(nperm_vector,gene_bin_theshold,genotype_df,sample_df,gene_list,bin_index,permutation_directory,TRUE,analysis_label, mc.cores = numCores)
    }
  }
}

#' returns the permutation based p values
return_permutation_p_vals_from_permutation_by_bin <- function(permutation_directory,analysis_label,nperm,gene_list, genotype_df, sample_df, numCores){
  permutation_p_numerator <- numeric(length = length(gene_list))
  estimate_true_vector <- numeric(length = length(gene_list))
  not_permuted_list <- list()
  
  for(bin_index in 1:length(gene_list)){
    nperm_vector <- 1:nperm
    print(paste0("Bin Index ", bin_index, " out of ", length(gene_list)))
    success <- NULL
    files <- paste0(permutation_directory,analysis_label,"_bin_", bin_index, "_perm_", 1:nperm,"_stat.RDS")
    # success <- as.numeric(gsub(paste0(analysis_label, "_bin_", bin_index,"_perm_|_stat.RDS"), "", files))
    
    # files <- paste0(permutation_directory, analysis_label_var, "_perm_", 1:nperm,"_stat.RDS")
    print("Loading in files...")
    temp_stat <- mclapply(files, readRDS, mc.cores = numCores)
    print("Extracting P values...")
    p_vals_of_random <- sapply(1:nperm, function(x) temp_stat[[x]][[1]]$p.value)
    
    not_permuted <- gene_bin_theshold(1, genotype_df, sample_df,gene_list,bin_index,permutation_directory,FALSE,analysis_label)
    permutation_p_numerator[bin_index] <- sum(not_permuted[[1]]$p.value >= p_vals_of_random)
    estimate_true_vector[bin_index] <- not_permuted[[1]]$estimate
    not_permuted_list[[length(not_permuted_list) + 1]] <- not_permuted
    
  }
  return_list <- list(permutation_p_numerator, permutation_p_numerator/nperm,estimate_true_vector,not_permuted_list)
  names(return_list) <- c("permutation_p_numerator","permutation_p_vals","estimate_true_vector","not_permuted_list")
  return(return_list)
}

#'' how many permutations do i need? method http://thompsonj.github.io/how-many-permutations
#' Function is icnomplete but can be expanded appropriately to provide guidance on the number of permutations necessary to achive a reasonable confidence interval
how_many_perms <- function(p, M,z){
  # debug data
  M <- 1000000
  z <- 1.645 #1.96 
  p <- 2.68 * 10^-5
  s = sqrt((p*(1-p)) / M)
  (e_min = (p-(z*s)))
  (e_max = (p+(z*s)))
  
  p <- 2.68 * 10^-5
  alpha <- 10^-5
  z <- 1.96
  epsilon <- 5*10^-5
  (M_needed <- (p*(1-p)) / (((alpha - epsilon - p)/z)^2))
  
  # M_needed <- (p*(1-p)) / (((e_max - p)/z)^2)
}


# return_stat_for_bin_size <- function(dummy_var,size_of_bin, genotype_df, sample_df, permutation_directory, analysis_var, permutation_var){
#   genotype_unique_sorted_df <- (genotype_df[order(genotype_df$gnomAD.Gene.oe_lof_upper),] %>% 
#                                   distinct(Sample.Name, Gene.Name, .keep_all = TRUE) %>% 
#                                   distinct(Sample.Name, .keep_all = TRUE)) %>% 
#     select(Sample.Name, Gene.Name, epi_phenotype, Sample.Phenotype,gnomAD.Gene.oe_lof_upper)
#   print(size_of_bin)
#   end_index <- nrow(genotype_unique_sorted_df) - size_of_bin+1
#   gene_list <- lapply(1:end_index,function(x) genotype_unique_sorted_df$Gene.Name[x:(x+size_of_bin-1)])
#   one_bin_size_list <- lapply(1:length(gene_list), 
#                               function (x) gene_bin_theshold(dummy_var, genotype_df, sample_df,gene_list,
#                                                              x,permutation_directory, paste0(analysis_var,"_", x), analysis_var))
#   return(one_bin_size_list)
# }
# 
# 
#CDF functions-----------

#' returns a df with info about % excess
#' return_total_samples<-function(atav_filtered_df,sample_df){
#'   total_samples <- sample_df %>% group_by(epi_phenotype) %>% summarize(total_subj_var = n())
#'   total_samples <- total_samples %>% rowwise() %>% mutate(variant_var=nrow(atav_filtered_df[atav_filtered_df$epi_phenotype==epi_phenotype,] %>% distinct(Sample.Name, .keep_all=TRUE)))
#'   total_samples <- total_samples %>% mutate(qv_percent = variant_var/total_subj_var)
#'   (total_samples <- total_samples %>% rowwise() %>% mutate(control_qv_fraction =  (total_samples %>% filter(epi_phenotype == "ctrl"))$qv_percent/qv_percent ))
#'   total_samples <- total_samples %>% rowwise() %>% mutate(benign_subjects = ceiling(total_subj_var * control_qv_fraction), benign_variants = ceiling(variant_var * control_qv_fraction))
#'   
#'   return(total_samples)
#' }
#' 
#' #' Returns CDF 
#' #' @param atav_temp
#' #' @param sample_temp
#' #' @param combo_values
#' #' @param total_samples
#' #' @param score_var
#' #' @param epi_prefix
#' 
#' return_cdf<-function(atav_temp, sample_temp, combo_values, total_samples){
#'   
#'   # generate cdfs
#'   combo_values_df <- data.frame("score"=combo_values)
#'   pair1_total <- (total_samples %>% filter(epi_phenotype == epi_prefix[1]))$total_subj_var
#'   pair2_total <- (total_samples %>% filter(epi_phenotype == epi_prefix[2]))$total_subj_var
#'   ctrl_total <- (total_samples %>% filter(epi_phenotype == "ctrl"))$total_subj_var
#'   
#'   
#'   ctrl_ecdf_obs<-ecdf(unlist((atav_temp %>% filter(epi_phenotype=="ctrl"))[score_var]))
#'   pair1_ecdf_obs<-ecdf(unlist((atav_temp %>% filter(epi_phenotype==epi_prefix[1]))[score_var]))
#'   pair2_ecdf_obs<-ecdf(unlist((atav_temp %>% filter(epi_phenotype==epi_prefix[2]))[score_var]))
#'   
#'   
#'   # *********
#'   # this takes the weighted average but the weights are defined by the total weights
#'   # *********
#'   if (tolower(static_dynamic) == "static"){
#'     # print("Using static weights")
#'     pair1_ecdf_exp <- pair1_ecdf_obs(combo_values_df$score)*(1- (total_samples %>% filter(epi_phenotype == epi_prefix[1]))$control_qv_fraction) + ctrl_ecdf_obs(combo_values_df$score)*(total_samples %>% filter(epi_phenotype == epi_prefix[1]))$control_qv_fraction
#'     pair2_ecdf_exp <- pair2_ecdf_obs(combo_values_df$score)*(1- (total_samples %>% filter(epi_phenotype == epi_prefix[2]))$control_qv_fraction) + ctrl_ecdf_obs(combo_values_df$score)*(total_samples %>% filter(epi_phenotype == epi_prefix[2]))$control_qv_fraction
#'     
#'   } else if (tolower(static_dynamic) == "dynamic") {
#'     # print("Using dynamic weights")
#'     # *********
#'     # this takes the weighted average at that point
#'     # *********
#'     combo_values_df <- combo_values_df %>% 
#'       rowwise() %>% 
#'       mutate(pair1_sample_incl = nrow(atav_temp %>% filter(epi_phenotype==epi_prefix[1], eval(parse(text=score_var)) <= score) %>% select(Sample.Name) %>% distinct())/pair1_total, 
#'              pair2_sample_incl = nrow(atav_temp %>% filter(epi_phenotype==epi_prefix[2], eval(parse(text=score_var)) <= score) %>% select(Sample.Name) %>% distinct())/pair2_total,
#'              ctrl_sample_incl = nrow(atav_temp %>% filter(epi_phenotype=="ctrl", eval(parse(text=score_var)) <= score) %>% select(Sample.Name) %>% distinct())/ctrl_total)
#'     
#'     combo_values_df <- combo_values_df %>% 
#'       rowwise() %>% 
#'       mutate(pair1_sample_phi = 1- (ctrl_sample_incl/(pair1_sample_incl )),
#'              pair2_sample_phi = 1- (ctrl_sample_incl/(pair2_sample_incl)))
#'     
#'     pair1_ecdf_exp <- pair1_ecdf_obs(combo_values_df$score)*combo_values_df$pair1_sample_phi + ctrl_ecdf_obs(combo_values_df$score)*(1-combo_values_df$pair1_sample_phi)
#'     pair2_ecdf_exp <- pair2_ecdf_obs(combo_values_df$score)*combo_values_df$pair2_sample_phi + ctrl_ecdf_obs(combo_values_df$score)*(1-combo_values_df$pair2_sample_phi)
#'   } else {
#'     error("static_dynamic has a bad value")
#'   }
#'   
#'   
#'   return(display_df <-data.frame("x_vals"= rep(combo_values,6), "x_vals_scaled" = rep(combo_values,6)/max(combo_values), "ecdf_val" = c(ctrl_ecdf_obs(combo_values),ctrl_ecdf_obs(combo_values), pair1_ecdf_obs(combo_values), pair1_ecdf_exp, pair2_ecdf_obs(combo_values), pair2_ecdf_exp), "exp_obs"=rep(c(rep("Observed",length(combo_values)),rep("Excess",length(combo_values))),3), "epi_phenotype"=rep(c("ctrl",epi_prefix),each=(length(combo_values)*2))))
#'   
#' }
#' 
