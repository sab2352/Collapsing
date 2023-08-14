# fp_create_gene_set.R: creates a new collapsing model using a gene set to subset an existing model

'Usage:
  fp_create_gene_set.R --model=<model> --resolution=<resolution> --cluster=<cluster> --gene_set_path=<gene_set_path> --gene_set_name=<gene_set_name> --ccds_genes_path=<ccds_genes_pathj> [--debug=<debug>]

  Options:
  -h --help
    --model=<model> model to extract variants using gene set
    --resolution=<resolution> resolution used to create clusters
    --cluster=<cluster> cluster to analyze
    --gene_set_path=<gene_set_path> path to gene set file. should be a .txt file. no header. First line of gene_set text file should start with # and be a description of the gene set. This line will be skipped when reading in file.
    --gene_set_name=<gene_set_name> name of gene set. no spaces.
    --ccds_genes_path=<ccds_genes_path> path to ccds genes in the atav database. will reject gene list and error if some genes are not in database
    --debug=<debug> [default: FALSE]
' -> doc

library(docopt)
library(here)
library(tidyverse)
library(data.table)
library(logr)
# library(lintr)
# library(styler)
# lintr::use_lintr(type = "tidyverse")

"%!in%" <- Negate("%in%")

arguments <- docopt(doc, version = "fp_create_gene_set.R 1.1")

if (arguments$debug == "TRUE") {
  arguments <- list()
  arguments$model <- "URFUNC_min80_igm"
  arguments$gene_set_path <- here("DefaultData", "gene_lists", "20180731_kidney_gu_genes.txt")
  arguments$ccds_genes_path <- here("DefaultData", "addjusted.CCDS.genes.index.r20.hg19.ensembl87.txt")
  arguments$gene_set_name <- "20180731_kidney_gu_genes"
  arguments$resolution <- "0_4"
  arguments$cluster <- "15"
}

#' loads in gene set and checks to make sure that the genes are in CCDS. Should probably also check that in ATAV but need to work on that.
#' @param arguments list of arguments that are passed from command line of program
#' @param list of genes to be checked
get_gene_set <- function(arguments, gene_set_pre_check) {
  ccds_genes_loaded_in <- fread(arguments$ccds_genes_path)
  ccds_genes_atav_format <- paste0("'", ccds_genes_loaded_in[[1]], "'")
  gene_set_genes <- paste0("'", gene_set_pre_check[[1]], "'")
  # "'ARSE'" %in% ccds_genes_atav_format
  # scratch
  # in_ccds_not_in_atav <- paste0("'",sort(ccds_genes_loaded_in$V1[ccds_genes_loaded_in$V1 %!in% atav_gene_symbols$V1]),"'")

  genes_in_gene_set_not_in_ccds <- sort(gene_set_genes[!(gene_set_genes %in% ccds_genes_atav_format)])

  if (length(genes_in_gene_set_not_in_ccds) > 0) {
    logr::log_print(genes_in_gene_set_not_in_ccds)
    logr::log_print("some genes included not in CCDS")
    return(NA)
  } else {
    logr::log_print("all genes found in CCDS")
    return(gene_set_genes)
  }
}

#' Initialize log file
#' @param arguments The arguments list with the scripts arguments
initialize_logfile <- function(arguments) {
  dir.create(here("Data"))
  dir.create(log_folder <- here("Data", "fp_create_gene_set_Log"))
  time_case_prefix <- paste0(gsub(":", "_", gsub(" ", "_", gsub("-", "_", Sys.time()))), "_")
  logr::log_open(paste0(log_folder, "/", time_case_prefix, "fp_create_gene_set_model_", arguments$model, "_cluster_", arguments$cluster, "_logfile.log"))
  logr::log_print(arguments)
}

initialize_logfile(arguments)
logr::log_print("log file initialized")

logr::log_print("Testing gene set to ensure compliance")
gene_set_pre_check <- fread(here(arguments$gene_set_path), header = FALSE, skip = 1)
gene_set <- get_gene_set(arguments, gene_set_pre_check)

if (is.na(gene_set[1])) {
  logr::log_print(gene_set)
  logr::log_print("Gene set not compliant")
  logr::log_print("error in gene set")
  stop()
  error()
  return()
  logr::log_print("shouldn't be here")
} else {
  logr::log_print("Gene set compliant")
  genotype_name <- list.files(path = here("Results", "Collapsing", paste0("LClust_res_", arguments$resolution, "_cluster_", arguments$cluster, "_FlashColl_07"), arguments$model), pattern = paste0(arguments$model, "_genotypes.csv$"), full.names = TRUE)
  nrow(genotype_df <- fread(genotype_name, header = TRUE))
  dir.create(save_path <- here("Results", "Collapsing", paste0("LClust_res_", arguments$resolution, "_cluster_", arguments$cluster, "_FlashColl_07/", arguments$model, "_", arguments$gene_set_name)))

  nrow(genotype_gene_set_df <- genotype_df %>% filter(`Gene Name` %in% gene_set))

  if (nrow(genotype_gene_set_df) == 0) {
    logr::log_print("There are 0 QVs in the loaded file")
    fwrite(x = data.frame(), file = here(save_path, paste0(arguments$model, "_", arguments$gene_set_name, "_matrix.txt.gz")), compress = "gzip", col.names = TRUE)
  }
  logr::log_print(names(genotype_gene_set_df))
  fwrite(x = genotype_gene_set_df, file = here(save_path, paste0(arguments$model, "_", arguments$gene_set_name, "_genotypes1.csv")), col.names = TRUE)
}
logr::log_print(names(genotype_gene_set_df))
fwrite(x = genotype_gene_set_df, file = paste0(save_path, "/",arguments$model,"_",arguments$gene_set_name,"_genotypes1.csv"), col.names = TRUE)

logr::log_print("finished")
logr::log_close()
