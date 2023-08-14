library(here)
library(docopt)
library(data.table)
library(logr)
library(tidyverse)
"%!in%" <- Negate("%in%")


nrow(gtex_df <- fread(here("DefaultData","GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")))
nrow(ccds_df <- fread(here("DefaultData","addjusted.CCDS.genes.index.r20.hg19.ensembl87.txt")))
nrow(atav_df <- fread(here("DefaultData","atavdb_unique_gene_072022.txt"), header = FALSE))

length(atav_ccds <- intersect(atav_df$V1, ccds_df$V1))
nrow(gtex_in_atav_ccds <- gtex_df %>% filter(Description %in% atav_ccds))
length(ccds_atav_genes_not_matched_in_gtex_df <- atav_ccds[atav_ccds %!in% gtex_df$Description])

nrow(bladder_genes <- gtex_in_atav_ccds %>% filter(Bladder > 50))
nrow(skin_genes <- gtex_in_atav_ccds %>% filter(`Skin - Not Sun Exposed (Suprapubic)` > 50))
nrow(Peripheral_nerve_genes <- gtex_in_atav_ccds %>% filter(`Nerve - Tibial` > 50))

write.table(file = here("Input","bladder_genes.txt"), x = bladder_genes %>% select(`# Bladder genes` = Description), quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(file = here("Input","skin_genes.txt"), x = skin_genes %>% select(`# Genes expressed more than 50 TPM in skin per GTEX` = Description), quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(file = here("Input","tibial_nerve_genes.txt"), x = Peripheral_nerve_genes %>% select(`# Genes expressed more than 50 TPM in tibial nerves per GTEX` = Description), quote = FALSE, row.names = FALSE, col.names = TRUE)
