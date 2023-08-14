library(tidyverse)
library(data.table)
library(broom)
library(openxlsx)

readGenos <- function(model, dir, nclust) {
  genos <- list()
  for (cl in nclust) {
    f <- list.files(paste0("Results/Collapsing", gsub("_", "", dir), "/LClust_res_", resolution, "_cluster_", cl, "_", dir, "FlashColl_07/", model), "genotypes.csv", full.names = TRUE)
    print(f)
    geno <- read_csv(f) %>% mutate(cl = cl)
    genos[[as.character(cl)]] <- geno
  }
  return(genos)
}

readGenosRecessive <- function(model, dir, nclust) {
  genos <- list()
  for (cl in nclust) {
    f <- list.files(paste0("Results/Collapsing", gsub("_", "", dir), "/LClust_res_", resolution, "_cluster_", cl, "_", dir, "FlashColl_07/", model), "comphet.csv", full.names = TRUE)
    print(f)
    geno <- fread(f, fill = TRUE, sep = ",")
    
    geno1 <- geno %>% select(contains("(#1)"))
    geno2 <- geno %>% select(contains("(#2)"))
    colnames(geno1) <- gsub(" \\(\\#1\\)", "", colnames(geno1))
    colnames(geno2) <- gsub(" \\(\\#2\\)", "", colnames(geno2))
    geno_new <- rbind(geno1, geno2) %>% distinct()
    geno_new <- geno_new %>% filter(!`Sample Name` == "") %>% mutate(cl = cl)
    genos[[as.character(cl)]] <- geno_new
  }
  return(genos)
}

getGeneGenos <- function(genos, gene) {
  
  gene <- paste0("'", gene, "'")
  geneGenosL <- lapply(genos, function (x) x %>% filter(`Gene Name` == gene))
  geneGenos <- do.call(rbind, geneGenosL)
                      
  return(geneGenos)
}

getGeneListGenos <- function(model, dir, nclust, gene.list, recessive = FALSE) {
  
  print(model)
  if (recessive) {
    genos <- readGenosRecessive(model, dir, nclust)
  } else {
    genos <- readGenos(model, dir, nclust)
  }
  
  geneGenosL <- lapply(gene.list, function (x) getGeneGenos(genos, x))

  return(geneGenosL)
}
