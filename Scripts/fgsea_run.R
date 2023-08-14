# fgsea_run.R: runs fgsea. Installation of fgsea library here https://bioconductor.org/packages/release/bioc/html/fgsea.html install reactome.db https://bioconductor.org/packages/release/data/annotation/html/reactome.db.html also need org.Hs.eg.db https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html

'Usage: 
  fgsea_run.R --resolution_var=<resolution_var> --model=<model> --min_sample=<min_sample> [--case_group=<case_group>] [--debug=<debug>] 
  
  Options:
  -h --help
  --resolution_var=<resolution_var> resolution for clustering.
  --model=<model> will run on a model if specified. If specified, no need for gene_list_case_path or gene_list_ctrl_path
  --min_sample=<min_sample> minimum size of case/control allowed to include cluster into combined cluster.
  --case_group=<case_group> [default: temp_group]
  --debug=<debug> [default: FALSE]

' -> doc

library(tidyverse)
library(data.table)
library(here)
library(docopt)
library(logr)
library(fgsea)
library(ggplot2)
library(readxl)
library(reactome.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

arguments <- docopt(doc, version = 'fgsea_run.R 1.0')
# debugging
if(arguments$debug == "TRUE"){
  arguments <- list()
  arguments$cluster <- "11"
  arguments$resolution_var <- "0_4"
  arguments$model <- "FlexFunc"
  arguments$min_sample <- "15"
  arguments$case_group <- "debug"
} 

time_case_prefix <- paste0(gsub(":", "_", gsub(" ", "_", gsub("-", "_", Sys.time()))), "_", arguments$case_group, "_",arguments$model, "_")

#' Initialize log file
#' @param arguments The arguments list with the scripts arguments
initialize_logfile <- function(arguments, label_var) {
  dir.create(here("Data"))
  dir.create(log_folder <- here("Data", paste0(label_var,"_Log")))
  logr::log_open(here(log_folder, paste0(time_case_prefix, label_var,"_logfile.log")))
  logr::log_print(arguments)
  logr::log_print(R.version)
}
# Create Logfile----
tryCatch({
  initialize_logfile(arguments, "fgsea_run")
}, error = function(e) {
  message("Caught an error creating log file: ", e$message)
})

# create rank file
tryCatch({
  dir.create(here("Results","fgsea"))  
  # adapt to your filename
  print(cl_size_path <- list.files(here("Results","Plots"), pattern = paste0("*lclustering_res_", arguments$resolution_var, "_cluster_sizes.txt$"), full.names = TRUE))
  cl_sizes <- fread(cl_size_path)
  nclust <- (cl_sizes %>% filter(case >= as.numeric(arguments$min_sample) & control >= as.numeric(arguments$min_sample)))$cluster %>% print()
  
  xlsx <- here("Results","CMH",paste0("CMH_exact_summary_lclust_res_", arguments$resolution_var, "_min_sample_", arguments$min_sample, ".xlsx"))
  cmhs <- xlsx %>%
    excel_sheets() %>%
    purrr::set_names() %>%
    map(read_excel, path = xlsx)
  nrow(cmh_model_case_only <- cmhs[[arguments$model]] %>% filter(is.na(estimate) | estimate > 1, p.value < 0.1)) #, p.value < 0.05
  # nrow(cmh_model_case_only <- cmhs[[arguments$model]] %>% filter(estimate < 1, p.value < 0.1)) #, p.value < 0.05
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'DYX1C1'"] <- "'DNAAF4'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'WDR63'"] <- "'DNAI3'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'WDR78'"] <- "'DNAI4'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'MAATS1'"] <- "'CFAP91'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'KIAA1024L'"] <- "'MINAR2'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'FAM84A'"] <- "'LRATD1'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'KIAA0355'"] <- "'GARRE1'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'C6orf211'"] <- "'ARMT1'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'FAM46C'"] <- "'TENT5C'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'TAZ'"] <- "'TAFAZZIN'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'MKL1'"] <- "'MRTFA'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'MFSD4'"] <- "'MFSD4A'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'C7orf72'"] <- "'SPATA48'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'HIST1H4E'"] <- "'H4C5'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'C2orf53'"] <- "'PRR30'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'VIMP'"] <- "'SELENOS'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'KIAA0226L'"] <- "'RUBCNL'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'PET112'"] <- "'GATB'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'B3GNT1'"] <- "'B3GNT2'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'HIST1H1A'"] <- "'H1-1'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'SRCRB4D'"] <- "'SSC4D'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'C2orf53'"] <- "'PRR30'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'IMPAD1'"] <- "'BPNT2'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'FAM71F2'"] <- "'GARIN1A'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'FAM65C'"] <- "'RIPOR3'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'FGFR1OP'"] <- "'FGFR1OP2'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'T'"] <- "'TBXT'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'EFCAB4B'"] <- "'CRACR2A'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'WDR65'"] <- "'CFAP57'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'HIST1H2AB'"] <- "'H2AC4'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'C2orf44'"] <- "'WDCP'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'AAED1'"] <- "'PRXL2C'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'EFCAB1'"] <- "'CLXN'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'MLK4'"] <- "'MAP3K21'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'C1orf228'"] <- "'ARMH1'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'C2orf54'"] <- "'MAB21L4'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'SEPT2'"] <- "'SEPTIN2'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'PRSS45'"] <- "'PRSS45P'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'SEPT4'"] <- "'SEPTIN4'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'FAM63A'"] <- "'MINDY1'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'ATP5A1'"] <- "'ATP5F1A'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'KIAA0556'"] <- "'KATNIP'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'ADCK3'"] <- "'COQ8A'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'KIAA1468'"] <- "'RELCH'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'KIAA1009'"] <- "'CEP162'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'ZUFSP'"] <- "'ZUP1'" 
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'INADL'"] <- "'PATJ'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'KIAA1524'"] <- "'CIP2A'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'WISP3'"] <- "'CCN6'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'TMEM2'"] <- "'CEMIP2'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'PDDC1'"] <- "'GATD1'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'KIAA1244'"] <- "'ARFGEF3'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'LEPREL1'"] <- "'P3H2'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'C1orf173'"] <- "'ERICH3'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'PCNXL2'"] <- "'PCNX2'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'LEPREL2'"] <- "'P3H3'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'FUK'"] <- "'FCSK'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'LRRC16B'"] <- "'CARMIL3'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'FBXO18'"] <- "'FBH1'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'CPSF3L'"] <- "'INTS11'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'LARS'"] <- "'LARS1'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'KIAA1279'"] <- "'KIFBP'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'IARS'"] <- "'IARS1'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'TTC37'"] <- "'SKIC3'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'GPR116'"] <- "'ADGRF5'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'C5orf42'"] <- "'CPLANE1'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'GCN1L1'"] <- "'GCN1'"
  cmh_model_case_only$gene[cmh_model_case_only$gene =="'MLLT4'"] <- "'AFDN'"
  # cmh_model_case_only$gene[cmh_model_case_only$gene =="'MRVI1'"] <- "'AFDN'"
  # cmh_model_case_only$gene[cmh_model_case_only$gene =="'HIST1H3G'"] <- "'AFDN'"
  # cmh_model_case_only$gene[cmh_model_case_only$gene =="'GPR113'"] <- "'AFDN'"
  # cmh_model_case_only$gene[cmh_model_case_only$gene =="'ADSSL1'"] <- "'AFDN'"
  # cmh_model_case_only$gene[cmh_model_case_only$gene =="'GPR133'"] <- "'AFDN'"
  # cmh_model_case_only$gene[cmh_model_case_only$gene =="'HIST3H3'"] <- "'AFDN'"
  # cmh_model_case_only$gene[cmh_model_case_only$gene =="'FAM154A'"] <- "'AFDN'"
  # cmh_model_case_only$gene[cmh_model_case_only$gene =="'RARS'"] <- "'AFDN'"
  # cmh_model_case_only$gene[cmh_model_case_only$gene =="'ASNA1'"] <- "'AFDN'"
  # cmh_model_case_only$gene[cmh_model_case_only$gene =="'CXorf64'"] <- "'AFDN'"
  # cmh_model_case_only$gene[cmh_model_case_only$gene =="'C11orf49'"] <- "'AFDN'"
  # cmh_model_case_only$gene[cmh_model_case_only$gene =="'HDGFRP2'"] <- "'AFDN'"
  fgsea_rank_df <- data.frame(ID = mapIds(org.Hs.eg.db, gsub("'","",cmh_model_case_only$gene), 'ENTREZID', 'SYMBOL'),
                              t = (-log10(cmh_model_case_only$p.value )))
  if(sum(is.na(fgsea_rank_df$ID)) > 0){
    logr::log_print("The following genes did not get matched")
    logr::log_print(row.names(fgsea_rank_df)[is.na(fgsea_rank_df$ID)])
    error()
    stop()
    return()
  }
  
  write.table(fgsea_rank_df,file=rnk.file <- here("Results","fgsea",paste0(time_case_prefix,"expression.rnk")),quote=F,sep="\t",row.names=F)

}, error = function(e) {
  message("Caught an error creating log file: ", e$message)
})

tryCatch({
  ranks <- read.table(rnk.file,
                      header=TRUE, colClasses = c("character", "numeric"))
  ranks <- setNames(ranks$t, ranks$ID)
  str(ranks)
  pathways <- reactomePathways(names(ranks))
  fgseaRes <- fgsea(pathways = pathways,
                    stats    = ranks,
                    minSize  = 15,
                    maxSize  = length(ranks) - 1)
  fgseaRes
  
  plot_to_save <- plotEnrichment(pathways[["Transport of small molecules"]],
                 ranks) + labs(title="Transport of small molecules")
  ggsave(plot = plot_to_save, filename = here("Results","fgsea",paste0(time_case_prefix, "plotEnrichment.pdf")), width = 3, height = 2, units = "in", dpi = 300)
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
                gseaParam=0.5)
  collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.05], 
                                        pathways, ranks)
  mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
    order(-NES), pathway]
  plot_to_save <- plotGseaTable(pathways[mainPathways], ranks, fgseaRes, 
                gseaParam = 0.5)
  ggsave(plot = plot_to_save, filename = here("Results","fgsea",paste0(time_case_prefix, "plot_fgsea_table.pdf")), width = 12, height = 1, units = "in", dpi = 300)
}, error = function(e) {
  message("Caught an error running fgsea: ", e$message)
})

tryCatch({
  save(fgseaRes, file = here("Results","fgsea",paste0(time_case_prefix,"fgsea_results.rds")))
  fwrite(fgseaRes, file=here("Results","fgsea",paste0(time_case_prefix,"fgsea_results.txt")), sep="\t", sep2=c("", " ", ""))
  fgseaResMain <- fgseaRes[match(mainPathways, pathway)]
  fgseaResMain[, leadingEdge := mapIdsList(
    x=org.Hs.eg.db, 
    keys=leadingEdge,
    keytype="ENTREZID", 
    column="SYMBOL")]
fwrite(fgseaResMain, file=here("Results","fgsea",paste0(time_case_prefix,"fgseaResMain.txt")), sep="\t", sep2=c("", " ", ""))
}, error = function(e) {
  message("Caught an error saving results: ", e$message)
})

tryCatch({
  logr::log_print("finished")
  logr::log_close()
}, error = function(e) {
  message("Caught an error closing log: ", e$message)
})
