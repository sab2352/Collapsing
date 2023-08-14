# lclust_Flash.R: creates clusters, Created by Gundula. Adapted for atav by Josh and others.

'Usage: 
  lclust_Flash.R [--case_group=<case_group>] [--resolution_var=<resolution_var>] [--debug=<debug>]
  
  Options:
  -h --help
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
  --resolution_var=<resolution_var> resolution for clustering. Typically 0.1-0.4 [default: 0.3]
  --debug=<debug> [default: FALSE]
  
' -> doc
library(tidyverse)
library(data.table)
library(here)
library(docopt)
library(cowplot)
library(gridExtra)
arguments <- docopt(doc, version = 'lclust_Flash.R 1.1')

if(arguments$debug == "TRUE"){
  arguments <- list()
  arguments$resolution_var <- "0.4"
  arguments$case_group <- "ic_20230501_singletons"
}
######### change this ############
dir <- ""

if (dir == "") {
    dir2 = ""
} else {
    dir2 <- paste0(dir, "_")
}



######### Initiate Log ############
dir.create(log_folder <- here("Data/lclust_FlashLog"))
time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_",arguments$case_group, "_")
logr::log_open(paste0(log_folder,"/", time_case_prefix,"lclust_Flash_logfile.log"))
logr::log_print(arguments)

logr::log_print(plot_path <- here("Results","Plots"))

if(!dir.exists(plot_path)) {
    dir.create(plot_path)
} else {
  logr::log_print(paste0("Directory ", plot_path, " already exists"))
}

# adapt paths to your structure
logr::log_print(dir)
logr::log_print(pc_path <- here(paste0("Results/KinshipFlashPCA/")))#used to have dir here but removed for now
logr::log_print(coverage_path <- here("Results/Coverage/"))
logr::log_print(evf <- list.files(pc_path, "*_flashpca_eigenvectors$", full.names = TRUE)[1])

pcs <- fread(evf)
pcs <- pcs %>% rename_at(vars(U1:U10), ~ paste0("PC", 1:10))

sf <- list.files(pc_path, "*_kinship_pruned_sample.txt$", full.names = TRUE)[1]
samples <- fread(sf, header = FALSE)

#### optional: remove low coverage samples ####
cf <- list.files(coverage_path, pattern = "_sample.summary.csv", full.names = TRUE)[1]
cov <- fread(cf, header = TRUE)
low_cov <- (cov %>% filter(`%Overall_Bases_Covered` < 0.85))$Sample

pcs <- pcs %>% filter(!IID %in% low_cov)
samples <- samples %>% filter(!V2 %in% low_cov)
####

cases <- samples %>% filter(V6 == 2)
controls <- samples %>% filter(V6 == 1)

##### include Ancestry from SequenceDB for plotting ####

######### change this ############
# file from Sequence with predicted ancestry information (either 1 file for all samples or split up into multiple files)
case_data <- list.files(here("Data/"), "*case.seq.csv$", full.names = TRUE) #
ctrl_data <- list.files(here("Data/"), "*ctrl.seq.csv$", full.names = TRUE)#
data_1 <- fread(case_data)# %>% filter(!(tolower(sample_internal_name) %in% tolower(exclude_list)))
data_2 <- fread(ctrl_data)
full_data <-rbind(data_1,data_2)
# data <- fread(paste0(sample_path, "sample_result_16940084_1561750729.csv"))
# 
# case_data <- fread(paste0(sample_path, "sample_result_16940084_1561752175.csv"))
# 
# full_data <- rbind(data, case_data)
# full_data <- data # option if all are in a single file (comment out case_data lines)
if ("CHGVID" %in% colnames(full_data)) {
    full_data <- full_data %>% rename(sample_internal_name = CHGVID) # if you have an old version of the file
}
full_data <- full_data[!duplicated(full_data$sample_internal_name),]

full_data <- full_data %>% mutate("Ancestry" = case_when(Caucasian_prob > 0.75 ~ "European",
                                            MiddleEastern_prob > 0.75 ~ "MiddleEastern",
                                            Hispanic_prob > 0.75 ~ "Latino",
                                            EastAsian_prob > 0.75 ~ "EastAsian",
                                            SouthAsian_prob > 0.75 ~ "SouthAsian",
                                            African_prob > 0.75 ~ "African",
                                            TRUE ~ "Admixed")) %>% 
    select(sample_internal_name, Ancestry, exomeKit)

pcs <- pcs %>% left_join(full_data, by = c("IID" = "sample_internal_name"))
pcs <- pcs %>% left_join(samples, by = c("IID" = "V2")) %>% 
    mutate(Phenotype = case_when(IID %in% cases$V2 ~ "case",
                                 IID %in% controls$V2 ~ "control",
                                 TRUE ~ "other"),
           Gender = case_when(V5 == 1 ~ "male",
                              V5 == 2 ~ "female",
                              TRUE ~ "other")) %>% 
    select(-c(V1:V8))

#### Louvain Clustering ####

library(Seurat)
library(Matrix)

total <- pcs %>% select(PC1:PC6) %>% as.matrix()
rownames(total) <- pcs$IID

total_nn <- FindNeighbors(total)
total_snn <- Matrix(total_nn[[2]], sparse=TRUE)

######### change this ############
resolution <- as.numeric(arguments$resolution_var) # 0.1-0.4
resolution2 <- gsub("\\.", "_", as.character(resolution))

total_clust <- FindClusters(total_snn, resolution = resolution)

saveRDS(total_clust, here(plot_path, paste0(dir2,time_case_prefix, "lclust_res_", resolution2, ".RDS")))

# total_clust <- readRDS(paste0(plot_path, "/", dir2, "lclust_res_", resolution2, ".RDS"))

pcs <- pcs %>% mutate(cluster = total_clust[,1])

#### visualize with UMAP ####

library(umap)

custom.settings <- umap.defaults
custom.settings$min_dist <- 0.1
custom.settings$random_state <- 188 # try another seed if the plot is not nice
custom.settings$verbose <- TRUE
custom.settings$metric <- "pearson"

umap <- umap(total, config = custom.settings)
saveRDS(umap, here(plot_path, paste0(dir2, "umap.RDS")))

y <- data.frame(umap$layout)
y <- y %>% mutate(clusters = total_clust[,1])

library(RColorBrewer)
ncolors <- max(length(unique(pcs$cluster)), length(unique(pcs$Ancestry)))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(ncolors)

cluster_plot <- ggplot(y, aes(x = X1, y = X2, color = clusters)) +
    geom_point(size = 0.5, alpha = 0.6) +
    theme(legend.title=element_blank()) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme_minimal() + 
    scale_color_manual(values = mycolors) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave(filename = here(plot_path, paste0(dir2,time_case_prefix, "UMAP_", resolution2, ".png")), 
       device = "png", width = 6, height = 5)

Ancestry <- factor(pcs$Ancestry)
ancestry_plot <- ggplot(y, aes(x = X1, y = X2, color = Ancestry)) +
    geom_point(size = 0.5, alpha = 0.6) +
    theme(legend.title=element_blank()) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme_minimal()+ 
    scale_color_manual(values = mycolors) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave(filename = here(plot_path, paste0(dir2,time_case_prefix, "Ancestry_UMAP.png")), 
       device = "png", width = 6, height = 5)


cl_ca_co <- pcs %>% group_by(cluster, Phenotype) %>% tally() %>% spread(Phenotype, n)

disp_table <- tableGrob(cl_ca_co, theme = ttheme_minimal(), rows = NULL)
theme_set(theme_cowplot(font_family = "Arial"))
temp<- cowplot::plot_grid(plotlist = list(ancestry_plot,cluster_plot,disp_table), nrow = 1,  rel_widths = c(1,1,.5) ) #rel_heights = c(.5, .5, 1),
temp <- temp + theme(text = element_text(family = "Arial")) +
  annotate("text",label = c("A","B","C"),x = c(0.01,0.33,66),y=rep(.95,3), size = 5, fontface ="bold")  #+
  # annotate("text",label = c("Ancestry","Clusters"),
  #          x = c(0.2,0.5),y=rep(.96,2), size = 4, color = "black", fontface ="bold")
ggsave(plot = temp, here(plot_path,paste0(dir2,time_case_prefix,"resolution_",resolution2,"_Ancestry_plus_clusters_plus_table.pdf")), width = 10, height = 4, units = "in",dpi = 300)



#### table of cluster sizes ####
write.table(cl_ca_co, here(plot_path, paste0(dir2, time_case_prefix, "lclustering_res_", resolution2, "_cluster_sizes.txt")), sep = "\t", row.names = FALSE, quote = FALSE)

cl_ca_co_ratio <- cl_ca_co %>% mutate(ratio = formatC(case/control, digits = 2))
write.table(cl_ca_co_ratio, paste0(plot_path, "/", dir2, time_case_prefix, "lclustering_res_", resolution2, "_cluster_sizes_ratios.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


#### write out cluster sample files ####
for (cl in levels(pcs$cluster)) {
    ind <- pcs %>% filter(cluster == cl)
    samples_cl <- samples %>% filter(V2 %in% ind$IID)
    write.table(samples_cl, sep = "\t", quote=FALSE, row.names = FALSE, col.names = FALSE, file=paste0(pc_path, "/flashPCA_lclustering_res_", resolution2, "_cluster_", cl, "_sample.txt"))
}


#################### one plot #####################

cl_anc <- pcs %>% group_by(cluster, Ancestry, Phenotype) %>% tally() %>% spread(Phenotype, n) %>% group_by(cluster) %>% top_n(1, control) %>%
    group_by(Ancestry) %>% mutate(AncestryNew=paste0(Ancestry, row_number())) %>% as.data.frame() %>% select(cluster, AncestryNew) %>% print()

cl_ca_co_anc <- cl_ca_co %>% left_join(cl_anc) %>%
    rename("ancestry" = "AncestryNew") %>% select(cluster, ancestry, case, control) %>% print()

write.table(cl_ca_co_anc, paste0(plot_path, "/", dir2, time_case_prefix, "lclustering_res_", resolution2, "_clusters_ancestry_sizes.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

cl_ca_co_anc_ratio <- cl_ca_co_anc %>% mutate(ratio = formatC(case/control, digits = 2)) %>% print()
write.table(cl_ca_co_anc_ratio, paste0(plot_path, "/", dir2,time_case_prefix, "lclustering_res_", resolution2, "_clusters_ancestry_sizes_ratios.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


pcs <- pcs %>% left_join(cl_anc %>% select(cluster, AncestryNew))

Ancestry <- factor(pcs$AncestryNew)
ggplot(y, aes(x = X1, y = X2, color = Ancestry)) +
    geom_point(size = 0.5, alpha = 0.6) +
    theme(legend.title=element_blank()) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme_minimal() +
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          legend.text = element_text(size=14),
          legend.title = element_text(size=14)) +
    scale_color_manual(values = mycolors) +
    guides(colour = guide_legend(override.aes = list(size = 6)))
ggsave(filename = paste0(plot_path, "/", dir2,time_case_prefix, "UMAP_", resolution2, "_Ancestry_comb.png"),
       device = "png", width = 8, height = 5)



#### group by gender ####

cl_ca_co_gender <- pcs %>% group_by(cluster, Phenotype, Gender) %>% tally() %>% spread(Phenotype, n) %>% as.data.frame()
write.table(cl_ca_co_gender, paste0(plot_path, "/", dir2,time_case_prefix, "lclustering_res_", resolution2, "_cluster_sizes_gender.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
cl_ca_co_gender_ratio <- cl_ca_co_gender %>% mutate(ratio = formatC(case/control, digits = 2))
write.table(cl_ca_co_gender_ratio, paste0(plot_path, "/", dir2,time_case_prefix, "lclustering_res_", resolution2, "_cluster_sizes_gender_ratios.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

#### write out cluster sample gender files ####
for (cl in levels(pcs$cluster)) {
    ind <- pcs %>% filter(cluster == cl)
    samplesCl <- samples %>% filter(V2 %in% ind$IID) %>% filter(V5 == 1)
    write.table(samplesCl, sep = "\t", quote=FALSE, row.names = FALSE, col.names = FALSE,
                file=paste0(pc_path, "/flashPCA_lclustering_res_", resolution2, "_cluster_", cl, "_sample_male.txt"))
    samplesCl <- samples %>% filter(V2 %in% ind$IID) %>% filter(V5 == 2)
    write.table(samplesCl, sep = "\t", quote=FALSE, row.names = FALSE, col.names = FALSE,
                file=paste0(pc_path, "/flashPCA_lclustering_res_", resolution2, "_cluster_", cl, "_sample_female.txt"))
}

logr::log_print("Finished")
logr::log_close()

return()
stop()
#haven't debugged past here
##################### SubGroups ################################


pcsAll <- pcs

dir <- "Adult"

samples_new <- fread(here(paste0("Data/2021-04-23CasesAndControls", dir, ".txt")))

if (dir == "") {
    dir2 = ""
} else {
    dir2 <- paste0(dir, "_")
}

plot_path <- here(paste0("Plots_", dir, "/"))
plot_path <- gsub("_/", "/", plot_path)

if(!dir.exists(plot_path)) {
    dir.create(plot_path)
} else {
    print((paste0("Directory ", plot_path, " already exists")))
}

pc_path <- here(paste0("Results/KinshipFlashPCA", dir))

if(!dir.exists(pc_path)) {
    dir.create(pc_path)
} else {
    print((paste0("Directory ", pc_path, " already exists")))
}


pcs <- pcs %>% filter(IID %in% samples_new$V2)
samples_pruned <- samples_new %>% filter(V2 %in% pcs$IID)
write.table(samples_pruned, sep = "\t", quote=FALSE, row.names = FALSE, col.names = FALSE,
            file=paste0(pc_path, "/", dir, "_kinship_pruned_sample.txt"))


#### table of cluster sizes ####
cl_ca_co <- pcs %>% group_by(cluster, Phenotype) %>% tally() %>% spread(Phenotype, n) %>% print()
write.table(cl_ca_co, paste0(plot_path, "/", dir2, "lclustering_res_", resolution2, "_cluster_sizes.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

cl_ca_co_ratio <- cl_ca_co %>% mutate(ratio = formatC(case/control, digits = 2))
write.table(cl_ca_co_ratio, paste0(plot_path, "/", dir2, "lclustering_res_", resolution2, "_cluster_sizes_ratios.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


#### write out cluster sample files ####
for (cl in levels(pcs$cluster)) {
    ind <- pcs %>% filter(cluster == cl)
    samples_cl <- samples %>% filter(V2 %in% ind$IID)
    write.table(samples_cl, sep = "\t", quote=FALSE, row.names = FALSE, col.names = FALSE,
                file=paste0(pc_path, "/flashPCA_lclustering_res_", resolution2, "_cluster_", cl, "_sample.txt"))

}

#### group by gender ####

cl_ca_co_gender <- pcs %>% group_by(cluster, Phenotype, Gender) %>% tally() %>% spread(Phenotype, n) %>% as.data.frame()
write.table(cl_ca_co_gender, paste0(plot_path, "/", dir2, "lclustering_res_", resolution2, "_cluster_sizes_gender.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
cl_ca_co_gender_ratio <- cl_ca_co_gender %>% mutate(ratio = formatC(case/control, digits = 2))
write.table(cl_ca_co_gender_ratio, paste0(plot_path, "/", dir2, "lclustering_res_", resolution2, "_cluster_sizes_gender_ratios.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

#### write out cluster sample gender files ####
for (cl in levels(pcs$cluster)) {
    ind <- pcs %>% filter(cluster == cl)
    samplesCl <- samples %>% filter(V2 %in% ind$IID) %>% filter(V5 == 1)
    write.table(samplesCl, sep = "\t", quote=FALSE, row.names = FALSE, col.names = FALSE,
                file=paste0(pc_path, "/flashPCA_lclustering_res_", resolution2, "_cluster_", cl, "_sample_male.txt"))
    samplesCl <- samples %>% filter(V2 %in% ind$IID) %>% filter(V5 == 2)
    write.table(samplesCl, sep = "\t", quote=FALSE, row.names = FALSE, col.names = FALSE,
                file=paste0(pc_path, "/flashPCA_lclustering_res_", resolution2, "_cluster_", cl, "_sample_female.txt"))
}

pcs <- pcsAll





###### additional plots (optional) ########
total_df <- total %>% as.data.frame()
Clusters <- as.factor(pcs$cluster)
Ancestry <- factor(pcs$Ancestry)

library(RColorBrewer)
ncolors <- max(length(unique(pcs$cluster)), length(unique(pcs$Ancestry)))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(ncolors)


plotPCs <- function(data, x, y, Clusters) {
    print(x)
    print(y)
    x1 <- rlang::sym(x)
    y1 <- rlang::sym(y)
    ggplot(data, aes(x = !!x1, y = !!y1, color = Clusters)) +
        geom_point(size = 0.5, alpha = 0.5) +
        theme(legend.title=element_blank()) +
        theme_minimal() +
        scale_color_manual(values = mycolors) +
        guides(colour = guide_legend(override.aes = list(size = 3)))
    ggsave(filename = paste0(plot_path, "/", dir2, "", x, "_", y, "_", resolution2, ".png"),
           device = "png", width = 8, height = 5)
}

plotPCsA <- function(data, x, y, Ancestry) {
    print(x)
    print(y)
    x1 <- rlang::sym(x)
    y1 <- rlang::sym(y)
    ggplot(data, aes(x = !!x1, y = !!y1, color = Ancestry)) +
        geom_point(size = 0.5, alpha = 0.5) +
        theme(legend.title=element_blank()) +
        theme_minimal() +
        scale_color_manual(values = mycolors) +
        guides(colour = guide_legend(override.aes = list(size = 3)))
    ggsave(filename = paste0(plot_path, "/", dir2, "Ancestry_", x, "_", y, ".png"),
           device = "png", width = 8, height = 5)
}

x <- paste0("PC", 1:6)
cb <- t(combn(x, 2))

map2(cb[,1], cb[,2], function(x, y) plotPCs(total_df, x, y, Clusters))
map2(cb[,1], cb[,2], function(x, y) plotPCsA(total_df, x, y, Ancestry))


