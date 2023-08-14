library(tidyverse)  
library(data.table)
library(broom)
library(openxlsx)

# adapt to your data structure
readSummary <- function(model, dir, nclust) {
    gsum <- list()
    for (cl in nclust) {
        print(cl)
        # might need modification
        f <- list.files(here(paste0("Results/Collapsing", gsub("_", "", dir), "/LClust_res_", resolution, "_cluster_", cl, "_", dir, "FlashColl_07/", model)), "summary.csv$", full.names = TRUE)
        gsum[[as.character(cl)]] <- fread(f)
    }
    return(gsum)
}

useCMH <- function(gsum, gene) {
    cat(gene)
    conf <- lapply(gsum, function (x) x %>% filter(`Gene Name` == gene) %>% 
                       select(`Qualified Case`, `Unqualified Case`, `Qualified Ctrl`, `Unqualified Ctrl`))
    confa <- simplify2array(lapply(conf, function(x) array(unlist(x), dim = c(2,2))))
    
    cmh <- mantelhaen.test(confa, exact = TRUE)
    
    confp <- lapply(gsum, function (x) x %>% filter(`Gene Name` == gene) %>% 
                        select(`Qualified Case`, `Unqualified Case`, `Qualified Ctrl`, `Unqualified Ctrl`, `Fet P`))
    
    confb <- do.call(rbind, conf)
    confb <- bind_rows(colSums(confb)) %>% 
        mutate(`%QV+ Case` = `Qualified Case`/(`Qualified Case` + `Unqualified Case`)*100, 
               `%QV+ Ctrl` = `Qualified Ctrl`/(`Qualified Ctrl` + `Unqualified Ctrl`)*100)
    
    tcmh <- as_tibble(c(tidy(cmh), confb, unlist(confp)))
    
    return(tcmh)
}

getCMHres <- function(model, dir, nclust) {
    print(model)
    gsum <- readSummary(model, dir, nclust)
    gene.names <- lapply(gsum, function(x) x$`Gene Name`)
    gene.names.r <- Reduce(union, gene.names)
    
    gsum <- lapply(gsum, function(x) x %>% full_join(as.tibble(gene.names.r), by = c("Gene Name" = "value")) %>% 
                       mutate(`Qualified Case` = case_when(is.na(`Qualified Case`) ~ 0, TRUE ~ as.double(`Qualified Case`)),
                              `Unqualified Case` = case_when(is.na(`Unqualified Case`) ~ max(x$`Qualified Case` + x$`Unqualified Case`), TRUE ~ (`Unqualified Case`)),
                              `Qualified Ctrl` = case_when(is.na(`Qualified Ctrl`) ~ 0, TRUE ~ as.double(`Qualified Ctrl`)),
                              `Unqualified Ctrl` = case_when(is.na(`Unqualified Ctrl`) ~ max(x$`Qualified Ctrl` + x$`Unqualified Ctrl`), TRUE ~ (`Unqualified Ctrl`))
                       ))    
    
    cmhs <- lapply(gene.names.r, function(i) useCMH(gsum, i))
    
    cmhs.tidy <- unnest(tibble(cmhs), cols = c(cmhs)) %>% 
        mutate(gene = gene.names.r) %>% arrange(p.value) %>% 
        select(-c(method, alternative, statistic)) %>% 
        select(gene, p.value, everything())
    
    results_path <- here(paste0("Results/CMH_", gsub("_", "",  dir), "/"))
    results_path <- gsub("_/", "/", results_path)
    
    saveRDS(cmhs.tidy, paste0(results_path, dir, model, "_res_", resolution, "_cluster_", paste(nclust, collapse = "_"), "_exact_CMH.RDS"))
    return(cmhs.tidy)
}
