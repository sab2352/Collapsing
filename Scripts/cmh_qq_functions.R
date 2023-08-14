library(tidyverse)
# library(dplyr)
# library(tidyr)
# library(tibble)
# library(readr)
# library(ggplot2)
# library(purrr)

library(broom)
library(readxl)
library(data.table)
library(ggrepel)
library(gridExtra)
library(cowplot)


readPermTable <- function(dir, model, perm) {
    print(perm)

    perm_file <- paste0(permut_path, dir, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_", perm, ".RDS")

     if(!file.exists(perm_file)) {
        if(!str_detect(model, "recessive|dominant")) {
            model <- paste0("dominant", model)
            perm_file <- paste0(permut_path, dir, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_", perm, ".RDS")
        }
        if(!file.exists(perm_file)) {
            perm_file <- paste0(permut_path, dir, model, "_CMH_perm_res_", resolution, "_min_sample_", min_sample, "_perm_", perm, ".RDS")
            if(!file.exists(perm_file)) {
            perm_file <- paste0(permut_path, dir, model, "_CMH_perm_res_", resolution, "_genderstrat_min_sample_", min_sample, "_perm_", perm, ".RDS")            
            }
        }
    }

    if(file.exists(perm_file)) {
        cmhp <- readRDS(perm_file)
        return(cmhp$p.value)
    } else {
        print(perm_file)
        return(NULL)
    }

}

#' plotQQ function displays and prints the QQ plots
plotQQ <- function(dir, model, cmh, nperm, top_plot, sig_line,  text_label_color_var, y_lim_low, title_size, case_group, label_size, mult_comp_method) {
    # create output folder
    dir.create(plot_dir <- here("Results","CMH",model))

    permRes <- lapply(nperm, function(x) readPermTable(dir, model, x))
    permRes.matrix <- do.call(cbind, permRes)
    permRes.matrix <- permRes.matrix[1:nrow(cmh),]
    pVal <- rowMeans(permRes.matrix)
    quant <- apply(permRes.matrix, 1, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
    
    obs <- cmh$p.value
    exp <- pVal

    obs <- obs[which(!is.na(obs))]
    exp <- exp[which(!is.na(exp))]
    obs[obs > 1] <- 1
    exp[exp > 1] <- 1
    exp <- exp[1:length(obs)]
    
    gws <- 0.05 / length(obs)
    reg_obs <- obs[(obs > gws) & (exp > gws) & (exp < 1) & (obs < 1)]
    reg_exp <- exp[(obs > gws) & (exp > gws) & (exp < 1) & (obs < 1)]
    reg_obs <- qchisq(reg_obs, 1, lower.tail = FALSE)
    reg_exp <- qchisq(reg_exp, 1, lower.tail = FALSE)
    
    obs <- -log10(obs)
    exp <- -log10(exp)
    
    uPerc <- -log10(quant[1,])[1:length(obs)]
    lPerc <- -log10(quant[2,])[1:length(obs)]
    
    d <- data.frame(gsub("'", "", cmh$gene)[!is.na(cmh$p.value)], obs, exp, uPerc, lPerc)
    names(d) <- c("gene", "obs", "exp", "upper", "lower")
    
    lin.reg <- lm(reg_obs ~ 0 + reg_exp[1:length(reg_obs)])
    lambda <- lin.reg$coefficients
    
    ggplot(d, aes(exp, obs)) +
        geom_point(colour="red") +
        geom_text_repel(data = subset(d, obs > 6), aes(label = gene)) +
        geom_line(aes(exp, upper), colour="orange") +
        geom_line(aes(exp, lower), colour="green") +
        geom_abline(slope = 1, colour="blue") +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste0("QQ Plot: Observed vs. expected p-values. Lambda = ", round(lambda, digits = 4))) +
        xlab("Expected -log10(p)") +
        ylab("Observed -log10(p)")
    
    ggsave(filename = paste0(plot_dir, "/",case_group,"_",model, "_CMH_qq_exact_res_", resolution, "_min_sample_", min_sample, "_v1.png"), 
           device = "png", width = 8, height = 8)
    
    # josh style.
    library(cowplot)
    cmh$atav_gene <- cmh$gene
    # parameters for size/shape of plot. In the future, should be modifiable at command line
    theme_size <- 12
    labels_from_top <- 1.8
    gundi_format = TRUE
    axis_title_size <- 10

    case_enriched_genes <- gsub("'","",filter(cmh, is.na(estimate) | estimate > 1)$gene)
    text_label_color <- rep(text_label_color_var,10)

    d$atav_gene <- paste0("'",d$gene,"'")
    if(gundi_format){
        d$estimate <- cmh$estimate[which(d$atav_gene == cmh$atav_gene)]
    } else {
        d<-merge(d,cmh, by = "atav_gene", sort = FALSE)
        d$gene <- d$gene.x
        case_enriched_genes <- gsub("'","",filter(d, is.na(estimate) | estimate > 1)$gene)
    }
    d$case_ctrl <- ifelse(d$gene %in% case_enriched_genes,"Case Enriched","Control Enriched")
    point_color_range_fxn <- colorRampPalette(c("blue","yellow","red"))
    point_color_range <- point_color_range_fxn(7)
    d$dot_color_index <- case_when(is.na(d$estimate) | as.numeric(as.character(d$estimate)) >= 15 ~ 7,
                                   as.numeric(as.character(d$estimate)) >= 7.5 & as.numeric(as.character(d$estimate)) < 15 ~ 6,
                                   as.numeric(as.character(d$estimate))  > 1 & as.numeric(as.character(d$estimate)) < 7.5 ~ 5,
                                   as.numeric(as.character(d$estimate)) <= 1/15 ~ 1,
                                   as.numeric(as.character(d$estimate)) <= 1/7.5 & as.numeric(as.character(d$estimate)) > 1/15 ~ 2,
                                   as.numeric(as.character(d$estimate)) < 1 & as.numeric(as.character(d$estimate)) < 1/7.5 ~ 3,            
                                   TRUE ~ 4)
    d$dot_color <- point_color_range[d$dot_color_index]
    
    theme_set(theme_cowplot(font_family = "ArialMT"))
    if(max(select(d,upper,lower,exp,obs)) > top_plot) stop("The Y limit set by the user is too restrictive to accommodate the most significant p value. The upper y limit will need to be increased from value with --top_plot.")
    text_df <- data.frame(label_text = rev(subset(d, gene %in% case_enriched_genes[1:10])$gene),
                          label_x_vals = seq(from = 0.1, to= (max(d$exp) - 0.25), length.out = 10),
                          segment_end_x_start = rev(subset(d, gene %in% case_enriched_genes[1:10])$exp),
                          segment_end_y_start = rev(subset(d, gene %in% case_enriched_genes[1:10])$obs+0.1)) 
    text_df$label_text <- gsub("KIAA2022","NEXMIF",text_df$label_text)#This is done for epi25 paper
    plot_obj <- ggplot(d, aes(exp, obs)) +
        geom_point(aes(color=case_ctrl), color = d$dot_color, size = 0.5) +
        geom_line(aes(exp, upper), colour="green") +
        geom_line(aes(exp, lower), colour="green") +
        geom_abline(slope = 1, colour="black") +
        geom_text(data = text_df, aes(x = label_x_vals, y = rep(top_plot-labels_from_top,10), label=label_text), angle = 45, hjust = 0,size = label_size,color = rev(text_label_color)) +#
        geom_segment(data = text_df, aes(x = label_x_vals, xend = segment_end_x_start, y = rep(top_plot-labels_from_top-0.2,10), yend = segment_end_y_start), alpha = 0.25) +
        ggtitle(paste0("QQ Plot: Lambda = ", round(lambda, digits = 4))) +
        xlab(bquote('Expected: -Log'[10]~'(P-Value)')) +
        ylab(bquote('Observed: -Log'[10]~'(P-Value)')) +
        ylim(0,top_plot+.1) +
        theme_cowplot(font_size = theme_size) +
        theme(plot.title = element_text(hjust = 0.5, size = title_size), legend.position = "bottom", axis.title=element_text(size=axis_title_size), panel.background = element_rect(fill = "white", colour = "white"), plot.background = element_rect(fill = "white"))
    if(!is.na(sig_line))  plot_obj <- plot_obj +   geom_hline(yintercept = -1*log10(sig_line), linetype="dashed", color = "#CC79A7")
    
    ggsave(filename = paste0(plot_dir, "/",case_group,"_",model, "_CMH_qq_exact_res_", resolution, "_min_sample_", min_sample, "_v2.png"), 
           device = "png", width = 8, height = 8)
    
    #create txt file of genes for each models add command for each model to txt file
    write.table(subset(d, gene %in% case_enriched_genes[1:10])$gene, file = here(paste0("Results/CMH/", model, "_resolution_", resolution, "_min_sample_", arguments$min_sample, "_top_10_genes.txt")), row.names = FALSE, col.names = FALSE, quote = FALSE)
    cmmd <- paste("$atav --gene", paste0("Results/CMH/", model, "_resolution_", resolution, "_min_sample_", arguments$min_sample, "_top_10_genes.txt"), "--list-gene-anno --out", plot_dir)
    write_lines(cmmd, gene_anno, append = TRUE)

    # 3rd plot, this one has an attached table for top 10 case hits
    cmh$estimate[is.na(cmh$estimate)] <- Inf
    temp_table <- cmh %>% filter(gene %in% paste0("'", case_enriched_genes[1:10],"'")) %>% mutate(`Case w/wo` = sprintf("%i / %i",`Qualified Case`,`Unqualified Case`), `Ctrl w/wo` = sprintf("%i / %i",`Qualified Ctrl`,`Unqualified Ctrl`), `P Value` = sprintf("%.1e",p.value)) %>% select(Gene=gene, `P Value`, OR = estimate, `Case w/wo`, `Ctrl w/wo`)

    disp_table <- tableGrob(temp_table, theme = ttheme_minimal(6), rows = NULL)
    theme_set(theme_cowplot(font_family = "Arial"))
    temp<- cowplot::plot_grid(plot_obj,disp_table, nrow = 1, rel_widths = c(1,.6) ) #rel_heights = c(.5, .5, 1),,
    temp <- temp + theme(text = element_text(family = "Arial"))  #+
    ggsave(plot = temp, paste0(plot_dir, "/",case_group,"_",model, "_CMH_qq_exact_res_", resolution, "_min_sample_", min_sample, "_perm_", length(nperm) ,"_mc_",mult_comp_method,"_v3.pdf"), width = 8, height = 4, units = "in",dpi = 300)
    
    
    # if (mult_comp_method == "fdr") {
    #     obs <- p.adjust(cmh$p.value, method = "fdr", n = length(cmh$p.value)) 
    #     exp <- p.adjust(pVal, method = "fdr", n = length(pVal)) 
    # }
    # 
    # temp <- fread("/Volumes/jm4279/addjusted.CCDS.genes.index.r20.hg19.ensembl87.txt")
    # names(temp)[1] <- "gene"
    # names(temp)[2] <- "chr"
    # nrow(d)
    # nrow(temp)
    # nrow(temp_merge <- merge(d,temp, by= "gene", all.x = TRUE) %>% arrange(chr))
    # setdiff(temp_merge %>% select(gene,exp,obs,upper,lower), d)
    # setdiff( d, temp_merge %>% select(gene,exp,obs,upper,lower))
    # length(unique(d$gene))
    # length(unique(temp_merge$gene))
    # nrow(temp_merge)
    # duplicate(temp)
    # View(duplicate_df <- temp[duplicated(temp$gene),])
    # 
    # temp_merge$x_loc <- 1:nrow(temp_merge)
    # ggplot(temp_merge, aes(x=x_loc, y=obs)) + geom_point() + ylab(bquote('Observed: -Log'[10]~'(P-Value)')) +
        
    
}        

plotQQAll <- function(dir, model, cmhs, nperm, top_plot, sig_line, text_label_color, y_lim_low,title_size, case_group,label_size,mult_comp_method) {
    print(model)
    cmh <- cmhs[[model]]
    plotQQ(dir, model, cmh, nperm, top_plot, sig_line, text_label_color, y_lim_low,title_size, case_group,label_size, mult_comp_method)
}
