# library(tidyverse)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(here)

library(data.table)
library(broom)
library(parallel)

# adapt to your data structure
readMatrix <- function(model, dir, pc_path, cluster) {
    logr::log_print(sprintf("Start readMatrix model %s cluster %i", model,cluster))
    model_folder <-  here("Results","Collapsing",paste0("LClust_res_", resolution, "_cluster_", cluster, "_", dir, "FlashColl_07"),model)
    sample.file <- list.files(model_folder, "*_existing.sample.txt$",  full.names = TRUE)
    matrix.file <- list.files(model_folder, "*_matrix.txt.gz$", full.names = TRUE)
    
    ped <- fread(sample.file, header = FALSE, sep = "\t")
    samples <- ped$V2
    is.case <- as.matrix(ped[, 6] == 2)
    rownames(is.case) <- samples
    print(matrix.file)
    #    data <- fread(matrix.file, header = TRUE, sep = "\t")
    data <- read_delim(matrix.file, delim = "\t")
    data <- data[,-ncol(data)]
    data <- data %>% column_to_rownames("sample/gene") %>% as.matrix()
    data[data > 0] <- 1
    
    logr::log_print(sprintf("End readMatrix model %s cluster %i", model,cluster))
    return(list(data, is.case))
}

permuteMatrix <- function(matrix.data, perm) {
    # logr::log_print(sprintf("Start permute matrix %i", perm))
    data <- matrix.data[[1]]
    is.case <- matrix.data[[2]]
    n.samples <- length(is.case)
    n.cases <- length(which(is.case))
    n.ctrls <- n.samples - n.cases
    
    set.seed(perm)
    K <- sample.int(n.samples, size = n.cases, replace = FALSE)
    qv.cases <- rowSums(data[,K])
    qv.ctrls <- rowSums(data[,-K])
    
    gsum <- data.frame(qv.cases, n.cases - qv.cases, qv.ctrls, n.ctrls - qv.ctrls, stringsAsFactors = FALSE)
    gsum <- gsum %>% rownames_to_column(var = "Gene Name")
    colnames(gsum) <- c("Gene Name", "Qualified Case", "Unqualified Case", "Qualified Ctrl", "Unqualified Ctrl")
    # logr::log_print(sprintf("End permute matrix %i", perm))
    return(gsum)
}


useCMHp <- function(gsum, gene) {
#    print(gene)
    # logr::log_print(sprintf("start useCMHp for %s",gene))
    conf <- lapply(gsum, function (x) x %>% filter(`Gene Name` == gene) %>% 
                       select(`Qualified Case`, `Unqualified Case`, `Qualified Ctrl`, `Unqualified Ctrl`))
    confa <- simplify2array(lapply(conf, function(x) array(unlist(x), dim = c(2,2))))
#    print(confa)

    cmh <- mantelhaen.test(confa, exact = TRUE)
    tcmh <- as_tibble(c(tidy(cmh)))
    
    # logr::log_print(sprintf("end useCMHp for %s",gene))
    return(tcmh)
}

doPermutations <- function(matrix.data, model, perm) {
    # logr::log_print(sprintf("Start permute %s %i", model,perm))
    gsum <- lapply(matrix.data, function(x) permuteMatrix(x, perm))
    names(gsum) <- names(matrix.data)
    
    gene.names <- lapply(gsum, function(x) x$`Gene Name`)
    gene.names.r <- Reduce(union, gene.names)
    
    # logr::log_print(sprintf("For permute %s %i, starting gsum lapply", model,perm))
    gsum <- lapply(gsum, function(x) x %>% full_join(as.tibble(gene.names.r), by = c("Gene Name" = "value")) %>% 
                       mutate(`Qualified Case` = case_when(is.na(`Qualified Case`) ~ 0, TRUE ~ as.double(`Qualified Case`)),
                              `Unqualified Case` = case_when(is.na(`Unqualified Case`) ~ max(x$`Qualified Case` + x$`Unqualified Case`), TRUE ~ (`Unqualified Case`)),
                              `Qualified Ctrl` = case_when(is.na(`Qualified Ctrl`) ~ 0, TRUE ~ as.double(`Qualified Ctrl`)),
                              `Unqualified Ctrl` = case_when(is.na(`Unqualified Ctrl`) ~ max(x$`Qualified Ctrl` + x$`Unqualified Ctrl`), TRUE ~ (`Unqualified Ctrl`))
                       ))    
    
    # logr::log_print(sprintf("For permute %s %i, finish gsum lapply", model,perm))
    cmhs <- lapply(gene.names.r, function(i) useCMHp(gsum, i))
    # logr::log_print(sprintf("For permute %s %i, finish gene.names.r lapply", model,perm))
    
    cmhs.tidy <- unnest(tibble(cmhs)) %>% 
        mutate(gene = gene.names.r) %>% arrange(p.value) %>%  
        select(gene, p.value)
    
    # logr::log_print(sprintf("For permute %s %i, finish cmhs.tidy lapply", model,perm))
    permut_path <- here(paste0("Results/Permutations_", gsub("_", "",  dir), "/"))
    permut_path <- gsub("_/", "/", permut_path)
     
    logr::log_print(sprintf("save_path is %s",save_path <- paste0(permut_path, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_", perm, ".RDS")))

    saveRDS(object = cmhs.tidy, file = save_path)
    # logr::log_print(sprintf("End permute %s %i", model,perm))
}

getCMHresP <- function(model, dir, nclust, nperm = c(1), cores = 1) {
    
    logr::log_print(paste(model, nperm))
    logfile <- here("Data","cmh_permutation_lclust_Log",paste0("messagesperm_", dir, model, "_", nperm[1], "_", Sys.Date(), ".log"))
    con <- file(logfile, "w")
    sink(con, append = TRUE, type = "message")
    print(sprintf("In getCMHresP, analyzing model %s",model))
    
    matrix.data <- lapply(nclust, function(x) readMatrix(model, dir, pc_path, x))
    print("matrix.data loaded")
    names(matrix.data) <- as.character(nclust)
    print("matrix.data renamed")
    
    mclapply(nperm, function(x) doPermutations(matrix.data, model, x), mc.cores = cores)
    sink(type="message")
    close(con)
}
