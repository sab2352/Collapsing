# plot_umap_table.R: creates a figure showing the ancestry umap, cluster umap

'Usage: 
  plot_umap_table.R --resolution_var=<resolution_var> --min_sample=<min_sample> [--case_group=<case_group>]
  
  Options:
  -h --help
  --min_sample=<min_sample> minimum number of cases/controls of included clusters
  --resolution_var=<resolution_var> resolution for clustering.
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]

' -> doc

library(docopt)
library(ggpubr)
library(tidyverse)
library(data.table)
library(here)
library(logr)
library(RColorBrewer)
"%!in%" <- Negate("%in%")
arguments <- docopt(doc, version = 'plot_umap_table.R 1.1')


#' Returns a dataframe with the sample files and an additional column showing the cluster membership
cluster_sample_df <-function(resolution_var, nclust){
  logr::log_print(file_list <- here(paste0("Results/KinshipFlashPCA/flashPCA_lclustering_res_",resolution_var,"_cluster_",nclust,"_sample.txt")))
  read_add_cluster<-function(path,cluster){
    sample_df <- fread(path)
    sample_df$cluster <- cluster
    return(sample_df)
  }
  sample_list <- lapply(as.numeric(as.character(nclust)), function(x) read_add_cluster(file_list[x+1], x) )
  sample_combined_df <- do.call(rbind, sample_list)
  return(sample_combined_df)
  
}
ancestry_plus_cluster <- function(resolution_var, output_folder, seq_df_var, min_sample, plot_width = 12, plot_height = 6){
  
  # initiate paths and load in data
  plot_path <- here(paste0("Results/Plots/"))
  pc_path <- here(paste0("Results/KinshipFlashPCA"))
  evf <- list.files(pc_path, "*_flashpca_eigenvectors$", full.names = TRUE)[1]
  logr::log_print(nrow(pcs <- fread(evf)))
  pcs <- pcs %>% rename_at(vars(U1:U10), ~ paste0("PC", 1:10))
  # logr::log_print(nrow(sf <- list.files(pc_path, "*_kinship_pruned_sample.txt$", full.names = TRUE)))
  # logr::log_print(nrow(samples <- fread(sf, header = FALSE)))
  logr::log_print(cl_size_path <- list.files(plot_path, pattern = paste0("*lclustering_res_", resolution_var, "_cluster_sizes.txt$"), full.names = TRUE))
  cl_sizes <- fread(cl_size_path)
  logr::log_print(nclust <- (cl_sizes %>% filter(case >= min_sample & control >= min_sample))$cluster)
  logr::log_print(nrow(samples <- cluster_sample_df(resolution_var,cl_sizes$cluster)))
  
  names(seq_df_var)[names(seq_df_var) == "sample_internal_name"] <- "V1"
    
  logr::log_print(nrow(seq_sample_merge <- merge(samples, seq_df_var, by  = 'V1')))
    
  list_obj_list <- cluster_return_figure_list(seq_df_var, samples, pcs,resolution_var, min_sample)

  png(paste0(output_folder,"cluster_table.png"), width = 5000, height = 2100, res = 200)
  temp<-cowplot::plot_grid(plotlist =list_obj_list, labels = "AUTO",nrow=1)#align = "h",

  print(temp)
  dev.off()
  
  ggsave(plot = temp,filename = paste0(output_folder,"cluster_table.pdf"), width = plot_width, height = plot_height, units = "in", dpi = 400)
}

cluster_return_figure_list <- function(full_data, samples, pcs_orig,resolution_var, min_sample){
  cases <- samples %>% filter(V6 == 2)
  controls <- samples %>% filter(V6 == 1)

  plot_path <- here(paste0("Results/Plots/"))
  # full_data <- data # option if all are in a single file (comment out case_data lines)
  if ("V1" %in% colnames(full_data)) {
    full_data <- full_data %>% rename(sample_internal_name = V1) # if you have an old version of the file
  }
  full_data <- full_data[!duplicated(full_data$sample_internal_name),]
  perc_threshold <- 0.95
  full_data <- full_data %>% mutate("Geographic Ancestry" = case_when(Caucasian_prob > perc_threshold ~ "European",
                                                                      MiddleEastern_prob > perc_threshold ~ "Middle Eastern",
                                                                      Hispanic_prob > perc_threshold ~ "Latino",
                                                                      EastAsian_prob > perc_threshold ~ "East Asian",
                                                                      SouthAsian_prob > perc_threshold ~ "SouthAsian",
                                                                      African_prob > perc_threshold ~ "African",
                                                                      TRUE ~ "Admixed")) %>% 
    select(sample_internal_name, SubProject, `Geographic Ancestry`)
  
  pcs <- pcs_orig %>% left_join(full_data, by = c("IID" = "sample_internal_name"))
  pcs <- pcs %>% left_join(samples, by = c("IID" = "V2")) %>% 
    mutate(Phenotype = case_when(IID %in% cases$V2 ~ "case",
                                 IID %in% controls$V2 ~ "control",
                                 TRUE ~ "other"),
           Gender = case_when(V5 == 1 ~ "male",
                              V5 == 2 ~ "female",
                              TRUE ~ "other")) %>% 
    select(-c(V1:V8))
  
  # load in RDS data
  total_clust <- readRDS(list.files(plot_path, pattern = paste0("*lclust_res_",resolution_var,".RDS$"), full.names = TRUE))
  pcs <- pcs %>% filter(FID %in% rownames(total_clust))
  umap <- readRDS(here("Results/Plots/umap.RDS"))

  ncolors <- max(length(unique(pcs$cluster)), length(unique(pcs$Ancestry)))
  mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(ncolors)
  y <- data.frame(umap$layout)
  y <- y %>% mutate(Clusters = total_clust[,1])
  
  
  cluster_obj <- ggplot(y, aes(x = X1, y = X2, color = Clusters)) +
    geom_point(size = 0.5, alpha = 0.6) +
    theme(legend.title=element_blank()) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme_minimal() + 
    scale_color_manual(values = mycolors) +
    guides(colour = guide_legend(override.aes = list(size = 3)))+ theme(legend.text=element_text(size=15), axis.title = element_blank(),legend.title=element_text(size=15), axis.text = element_blank())
  # print(umap_obj)
  # ggsave(filename = paste0(plot_path, dir2, "UMAP_", resolution2, ".png"), 
  # device = "png", width = 6, height = 5)
  
  Ancestry <- factor(pcs$`Geographic Ancestry`)
  umap_obj <- ggplot(y, aes(x = X1, y = X2, color = Ancestry)) +
    geom_point(size = 0.5, alpha = 0.6) +
    theme(legend.title=element_blank()) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme_minimal()+ 
    # scale_color_manual(values = mycolors) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    theme(legend.text=element_text(size=15), axis.title = element_blank(),legend.title=element_text(size=15), axis.text = element_blank())
  # print(cluster_obj)
  # ggsave(filename = paste0(plot_path, dir2, "Ancestry_UMAP.png"), 
  # device = "png", width = 6, height = 5)

    (table_for_display <- samples %>% group_by(cluster, V6) %>% summarize(count_var = n()) %>% spread(V6, count_var, fill = 0) %>% select(Cluster = cluster, Cases = `2`, Controls = `1`))
    palette_var <- c(rep("pink",nrow(table_for_display)))
    palette_var[table_for_display$Cases >= min_sample & table_for_display$Controls >= min_sample] <- "white"
    
    t1 = ttheme(base_style = "default", 
                colnames.style = colnames_style(fill = "white",linecolor = "black", size = 17),
                tbody.style = tbody_style(fill = palette_var, linewidth = 1, linecolor = "black", size = 17)
    )
    p2<-ggtexttable(table_for_display, rows = NULL, theme = t1)#,padding = unit(c(3, 8), "mm")
    
  return(list(umap_obj,cluster_obj, p2))
} 

# ------main------
here(dir.create("Data"))
here(dir.create(output_folder <- "Results/Plots/ancestry with clusters/"))
time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_",arguments$case_group, "_")
logr::log_open(paste0(output_folder,"/", time_case_prefix,"plot_umap_table.log"))
logr::log_print(arguments)


dir.create(output_folder <- here(paste0("Results/Plots/ancestry with clusters/")))
case_seq_path <- list.files(here("Data"), pattern = "*case.seq.csv$", full.names = TRUE)
case_seq <- fread(case_seq_path)
ctrl_seq_path <- list.files(here("Data"), pattern = "*ctrl.seq.csv$", full.names = TRUE)
ctrl_seq <- fread(ctrl_seq_path)
ancestry_plus_cluster(arguments$resolution_var, output_folder, rbind(case_seq,ctrl_seq), as.numeric(arguments$min_sample))

# close log file
logr::log_close()
