# cmh_qq.R: Produces QQ plots. Created by Gundula. Adapted for atav by Josh and others.
# make sure you are in the right wd before loading here

'Usage: 
  cmh_qq.R --resolution_var=<resolution_var> --min_sample=<min_sample> --nperm_end=<nperm_end> [--top_plot=<top_plot>] [--sig_line=<sig_line>] [--text_label_color=<text_label_color>] [--y_lim_low=<y_lim_low>] [--title_size=<title_size>] [--case_group=<case_group>] [--label_size=<label_size>] [--mult_comp_method=<mult_comp_method>] [--debug=<debug>]
  
  Options:
  -h --help
  --min_sample=<min_sample> minimum number of cases/controls of included clusters
  --resolution_var=<resolution_var> resolution for clustering.
  --nperm_end=<nperm_end> number of permutations to incorporate 
  --top_plot=<top_plot> max y value displayed in qq plot [default: 8.9]
  --sig_line=<sig_line> p value to add horizontal line for significance [default: NA]
  --text_label_color=<text_label_color> color of text for label [default: black]
  --y_lim_low=<y_lim_low> i think this may never get used. will remove in future versions. [default: 6]
  --title_size=<title_size> Font size for title [default: 9]
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
  --label_size=<label_size> label_size is the font size of the gene labels for the top 10 hits [default: 3]
  --mult_comp_method=<mult_comp_method> [default: none]
  --debug=<debug> [default: FALSE]
  
  
' -> doc

library(docopt)
library(here)
source(here("Scripts/cmh_qq_functions.R"))
arguments <- docopt(doc, version = 'cmh_qq.R 1.1')

if(arguments$debug == "TRUE"){
  arguments<-list()
  arguments$min_sample <- "20"
  arguments$resolution_var <- "0_2"
  arguments$nperm_end <- "25"
  arguments$top_plot <- "8.9"
  arguments$sig_line <- "NA"
  arguments$text_label_color <- "black"
  arguments$y_lim_low <- "6"
  arguments$title_size <- "9"
  arguments$case_group <- "temp"
  arguments$mult_comp_method <- "fdr"
}

#########  ############
dir.create(log_folder <- here("Data","cmh_qq_Log"))
time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_")
logr::log_open(paste0(log_folder,"/", time_case_prefix,"cmh_qq_logfile.log"))
logr::log_print(arguments)


resolution <- arguments$resolution_var
min_sample <- as.numeric(arguments$min_sample)
nperm_end <- as.numeric(arguments$nperm_end)
sig_line <- ifelse(arguments$sig_line == "NA",NA,as.numeric(arguments$sig_line))
# text_label_color <- ifelse(arguments$text_label_color == "NULL",as.null("NULL"),arguments$text_label_color)


dir <- ""

if (dir == "") {
    dir = ""
} else {
    dir <- paste0(dir, "_")
}

plot_path <- paste0(here(paste0("Results/Plots_", gsub("_", "",  dir))), "/")
plot_path <- gsub("_/", "/", plot_path)
results_path <- paste0(here(paste0("Results/CMH_", gsub("_", "",  dir))), "/")
results_path <- gsub("_/", "/", results_path)
permut_path <- paste0(here(paste0("Results/Permutations_", gsub("_", "",  dir))), "/")
permut_path <- gsub("_/", "/", permut_path)


# adapt to your filename
print(cl_size_path <- list.files(plot_path, pattern = paste0("*lclustering_res_", resolution, "_cluster_sizes.txt$"), full.names = TRUE))
cl_sizes <- fread(cl_size_path)
nclust <- (cl_sizes %>% filter(case >= min_sample & control >= min_sample))$cluster %>% print()

xlsx <- paste0(results_path, dir, "CMH_exact_summary_lclust_res_", resolution, "_min_sample_", min_sample, ".xlsx")
# xlsx <- paste0(results_path, dir, "CMH_exact_summary_lclust_res_", resolution, "_min_sample_", min_sample, "_recessive.xlsx")
# xlsx <- paste0(results_path, dir, "CMH_exact_summary_lclust_res_", resolution, "_min_sample_", min_sample, "_genderstrat_recessive.xlsx")

cmhs <- xlsx %>%
    excel_sheets() %>%
    purrr::set_names() %>%
    map(read_excel, path = xlsx)

# tab names in xslx files can only be 30 characters
# names(cmhs)[grepl("dominantFlexiblePolyphen", names(cmhs))] <- "dominantFlexiblePolyphenDamaging"

#create file that contains list of atav commands
gene_anno <- paste0(results_path, "gene_anno_commands_resolution_", resolution, "_min_sample_", arguments$min_sample, ".txt")
file.create(gene_anno)

models <- names(cmhs) %>% print()

nperm <- 1:nperm_end

lapply(models, function(x) plotQQAll(dir, x, cmhs, nperm, top_plot = as.numeric(arguments$top_plot), sig_line=sig_line, text_label_color=arguments$text_label_color, y_lim_low=as.numeric(arguments$y_lim_low), title_size=as.numeric(arguments$title_size), case_group = arguments$case_group,
    label_size=as.numeric(arguments$label_size), mult_comp_method = arguments$mult_comp_method))


logr::log_print("Finished")
logr::log_close()

# model <- "dominantUltraRareEnsemble"
# cmh <- readRDS(paste0(dir, model, "_res_", resolution, "_exact_CMH.RDS"))
# plotQQ(dir, model, cmh, nperm)
