# fp_forest_plot.R: creates a forest plot, Created by Josh. Plotted on log scale. see discussion here. https://andrewpwheeler.com/2013/10/26/odds-ratios-need-to-be-graphed-on-log-scales/

'Usage:
  fp_forest_plot.R --output_directory=<output_directory> --path_to_fp_models=<path_to_fp_models> [--title_str=<title_str>]  [--pub_leg_pos=<pub_leg_pos>] [--pub_breaks=<pub_breaks>] [--pub_width=<pub_width>] [--pub_height=<pub_height>] [--x_val_factor=<x_val_factor>]  [--p_val_textsize=<p_val_textsize>]  [--p_line_space=<p_line_space>]     [--y_axis_alpha=<y_axis_alpha>]     [--palette_index=<palette_index>]  [--show_x_axis_title=<show_x_axis_title>]   [--y_top_buffer=<y_top_buffer>]     [--y_bottom_buffer=<y_bottom_buffer>]    [--p_alpha_p_val=<p_alpha_p_val>]     [--p_val_hjust=<p_val_hjust>]    [--pointrange_size=<pointrange_size>]     [--arrow_size=<arrow_size>]   [--sig_line_alpha=<sig_line_alpha>] [--show_per_w_p=<show_per_w_p>]     [--margin_vector=<margin_vector>]  [--y_max_pub=<y_max_pub>]  [--pub_width_table_mod=<pub_width_table_mod>] [--show_ci_w_p=<show_ci_w_p>] [--y_min_pub=<y_min_pub>] [--label_w_10_power=<label_w_10_power>]  [--show_or=<show_or>]  [--rel_width_vector=<rel_width_vector>] [--debug=<debug>] [--cores=<cores>] [--cowplot_font_size=<cowplot_font_size>]

  Options:
  -h --help
    --output_directory=<output_directory> path to output directory
    --path_to_fp_models=<path_to_fp_models> path to csv with models
    --title_str=<title_str> string for labelling output files [default: temp]
    --pub_leg_pos=<pub_leg_pos> Controls position of legend. Accepts none right left top bottom. also accepts vector for position within plot (e.g., 0.9,0.1). If using vector, no spaces after comma [default: none]
    --pub_breaks=<pub_breaks> Default labels on horizontal axis. Remember, log scale [default: .1,1,10]
    --pub_width=<pub_width> width in inches of pdf output [default: 8]
    --pub_height=<pub_height> height in inches of pdf output [default: 8]
    --x_val_factor=<x_val_factor> compresses or expands along the vertical axis [default: 1]
    --p_val_textsize=<p_val_textsize> font size for p value text [default: 3]
    --p_line_space=<p_line_space> the larger the number, the more space between the point and the p value label[default: 0.3]
    --y_axis_alpha=<y_axis_alpha> Sets alpha for vertical axes. Set to 0 and vertical axis line will disappear [default: 1]
    --palette_index=<palette_index>  index of color to be used for each line. must be a list of integers >=1 with , separation and no spaces [default: 1,2,3,4,5,6,7,8]
    --show_x_axis_title=<show_x_axis_title> [default: TRUE]
    --y_top_buffer=<y_top_buffer> Creates more space at top of plot. This is useful if a point has a p value label at the top and it is being cutoff. Increase site to allow more of top of plot to show [default: 0.5]
    --y_bottom_buffer=<y_bottom_buffer> [default: 0]
    --p_alpha_p_val=<p_alpha_p_val> value which will filter out point labels with OR and P. Set to -1 to see none, set to 2 to see all. Set to 0.05 and will only show labels for points with an uncorrected P < 0.05 [default: 0.05]
    --p_val_hjust=<p_val_hjust> Determines position of p value label relative to line. Range 0 - 1. Increase number to get p val text closer to vertical axis [default: 0.15]
    --pointrange_size=<pointrange_size> [default: 0.5]
    --arrow_size=<arrow_size> not currently used. likely delete [default: 0.5]
    --sig_line_alpha=<sig_line_alpha> [default: 0.2]
    --show_per_w_p=<show_per_w_p> If true, the percentage of cases/controls with a QV are shown with the p value annotation [default: TRUE]
    --margin_vector=<margin_vector> [default: 0.5,0,0,0]
    --y_max_pub=<y_max_pub> The maximum value on the horizontal axis. Its a y value because the coordinates are switched for display purposes [default: 10]
    --y_min_pub=<y_min_pub> The minimum value on the horizontal axis. Its a y value because the coordinates are switched for display purposes [default: 0.1]
    --pub_width_table_mod=<pub_width_table_mod> adjust table width relative to graph wide. e.g. 1 means will be one inch wider than forest plot [default: 0]
    --show_ci_w_p=<show_ci_w_p> show confidence intervals over points in forest plot [default: TRUE]
    --label_w_10_power=<label_w_10_power> if true, outputs x axis as 10^. If false, prints number. [default: TRUE]
    --show_or=<show_or> if true, shows OR label on top of data points [default: TRUE]
    --rel_width_vector=<rel_width_vector> relatives widths of final table with three panels [default: 1,.5,.1]
    --debug=<debug> if true, uses debugger arguments. only for developers [default: FALSE]
    --cores=<cores> numbers of cores for parallel rpocessing [default: 1]
    --cowplot_font_size=<cowplot_font_size> size of font in cowplot theme [default: 10]

' -> doc

library(docopt)
library(tidyverse)
library(data.table)
library(readxl)
library(here)
library(logr)
library(cowplot)
library(scales)
library(parallel)
source(here("Scripts", "fp_forest_plot_functions.R"))

arguments <- docopt(doc, version = "fp_forest_plot.R 1.1")

# debugging
if (arguments$debug == TRUE) {
  arguments <- list()
  # arguments$num_groups_var <- "Analysis"
  arguments$output_directory <- "/Volumes/projects/refractory_epilepsy/20221122_total/ClusteredCollapsing/Results/fp/"
  # arguments$left_border <- "-1"
  arguments$pub_width <- "4"
  arguments$pub_height <- "2"
  arguments$palette_index <- "1,2,3,4"
  arguments$group_color_palette_index <- "1,2,3,4"
  arguments$x_val_factor <- "1"
  # arguments$pub_textsize <- "5"
  arguments$pub_text_start_buf <- "0.1"
  arguments$label_height_factor <- "0.95"
  arguments$y_max_pub <- "11"
  # arguments$group_label_left_border <- "4"
  arguments$p_alpha_p_val <- "0.05"
  # arguments$color_margin <- "0.5"
  arguments$show_per_w_p <- "TRUE"
  arguments$p_line_space <- "0.3"
  arguments$p_val_hjust <- "0.15"
  arguments$p_val_textsize <- "3"
  arguments$show_x_axis_title <- "TRUE"
  arguments$pub_breaks <- "0.1,1,10"
  arguments$pointrange_size <- "0.5"
  arguments$sig_line_alpha <- "0.2"
  # arguments$dashed_alpha <- "0"
  arguments$y_top_buffer <- "0.5"
  arguments$y_bottom_buffer <- "0.5"
  # arguments$y_axis_width <- "0.8"
  arguments$y_axis_alpha <- "1"
  # arguments$x_axis_line_end_buffer <- "1"
  # arguments$x_axis_height <- "1.5"
  # arguments$group_text_size <- "3"
  arguments$group_label_alpha <- "1"
  arguments$margin_vector <- "0.5,0,0,0"
  arguments$pub_leg_pos <- "none"
  arguments$title_str <- "title_str2"
  arguments$log_flag <- "TRUE"
  # arguments$color_var <- "Analysis"
  # arguments$left_label_width <- ""
  arguments$arrow_size <- "0.5"
  arguments$show_ci_w_p <- "TRUE"
  arguments$group_label_vert_axis_alpha <- 1
  arguments$path_to_fp_models <- here("Input", "fp_models_refractory_nonrefractory_pgx.csv")
  arguments$cores <- "4"
  arguments$rel_width_vector <- "1,.5,.1"
  arguments$show_or <- "TRUE"
  arguments$label_w_10_power <- "TRUE"
  arguments$y_min_pub <- "0.1"
  arguments$cowplot_font_size <- "10"
}

time_case_prefix <- paste0(gsub(":", "_", gsub(" ", "_", gsub("-", "_", Sys.time()))), "_", arguments$title_str, "_")

#' creates a forest plot
#' @param arguments is a list with all of the command line arguments
display_fp <- function(arguments) {

  # Initialize variables
  palette_colors <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000")
  point_color_range_fxn <- colorRampPalette(c("blue", "yellow", "red"))
  point_color_range <- point_color_range_fxn(7)


  # convert variables
  logr::log_print(palette_index <- as.numeric(unlist(strsplit(arguments$palette_index, ","))))
  logr::log_print(pub_breaks <- as.numeric(unlist(strsplit(arguments$pub_breaks, ","))))
  logr::log_print(margin_vector <- as.numeric(unlist(strsplit(arguments$margin_vector, ","))))
  logr::log_print(rel_width_vector <- as.numeric(unlist(strsplit(arguments$rel_width_vector, ","))))
  p_line_space <- as.numeric(arguments$p_line_space)
  p_val_hjust <- as.numeric(arguments$p_val_hjust)
  p_val_textsize <- as.numeric(arguments$p_val_textsize)
  show_x_axis_title <- as.logical(arguments$show_x_axis_title)
  y_bottom_buffer <- as.numeric(arguments$y_bottom_buffer)
  y_top_buffer <- as.numeric(arguments$y_top_buffer)
  y_axis_alpha <- as.numeric(arguments$y_axis_alpha)
  num_groups_var <- "Analysis"
  output_directory <- arguments$output_directory
  arrow_size <- as.numeric(arguments$arrow_size)
  y_min_pub <- as.numeric(arguments$y_min_pub)
  log_flag <- as.logical(arguments$log_flag)
  label_w_10_power <- as.logical(arguments$label_w_10_power)

  color_var <- "Analysis"
  show_ci_w_p <- as.logical(arguments$show_ci_w_p)
  # *********
  # load in forest plot data
  # *********
  ifelse(!is.null(arguments[["fp_obj"]]), fisher_df <- arguments[["fp_obj"]], fisher_df <- fread(arguments$path_to_fp_data))
  logr::log_print(fisher_df)
  # *********


  if ((pub_width <- as.numeric(arguments$pub_width)) > 11) stop("pub width too large")
  if ((pub_height <- as.numeric(arguments$pub_height)) > 11) stop("pub height too large")
  logr::log_print(palette_colors <- palette_colors[palette_index])
  # *********
  # Correct for infinity and 0
  # *********
  fisher_df_original <- fisher_df
  non_inf_index <- fisher_df$OR != Inf & fisher_df$CI_high != Inf
  max_noninfinityor <- max(as.numeric(as.character(fisher_df$OR[non_inf_index])) + as.numeric(as.character(fisher_df$CI_high[non_inf_index])))
  if (length(fisher_df$CI_high[index <- fisher_df$OR == Inf]) > 0) fisher_df$CI_high[index] <- as.character(max_noninfinityor)
  if (length(fisher_df$OR[index <- fisher_df$OR == Inf]) > 0) fisher_df$OR[index] <- as.character(max_noninfinityor)

  # header and column info----
  fisher_df$x_val <- (1:nrow(fisher_df)) * (x_val_factor <- as.numeric(arguments$x_val_factor))

  header_str <- c("Condition/Gene", "Group", "W/ QV", "W/o QV", "% w/ QV", "OR", "P")
  text_column_start <- -16
  y_max <- 10
  condition_gene_width <- 5
  group_start <- text_column_start + condition_gene_width
  group_width <- 1.2
  w_qv_start <- group_start + group_width
  w_qv_width <- 1.2
  wo_qv_start <- w_qv_start + w_qv_width
  wo_qv_width <- w_qv_width
  percent_w_qv_start <- wo_qv_start + wo_qv_width
  percent_w_qv_width <- w_qv_width
  or_start <- percent_w_qv_start + percent_w_qv_width
  or_width <- w_qv_width
  p_start <- or_start + or_width
  p_width <- w_qv_width + 0.4

  y_axis_starts <- c(text_column_start, group_start, w_qv_start, wo_qv_start, percent_w_qv_start, or_start, p_start)
  line_space <- 0.45

  # Full figure, not publicatoin figure
  # ***********
  # Publication figure
  fisher_df_original <- fisher_df
  fisher_df$CI_low <- as.numeric(as.character(fisher_df$CI_low))
  fisher_df$CI_high <- as.numeric(as.character(fisher_df$CI_high))
  line_space <- 0
  y_max_pub <- as.numeric(arguments$y_max_pub)
  p_alpha <- ifelse(as.numeric(fisher_df$P_uncorrected) < (p_alpha_p_val <- as.numeric(arguments$p_alpha_p_val)), 1, 0)
  p_val_text <- ifelse(as.numeric(fisher_df$P_corrected) >= 0.01, sprintf("%.3f", as.numeric(fisher_df$P_corrected)), sprintf("%.1e", as.numeric(fisher_df$P_corrected)))
  p_val_text[p_val_text == "1.000"] <- "1"

  fisher_df$`Case % w/ QV` <- sprintf("%.1f", fisher_df$`Case w QV` / (fisher_df$`Case w QV` + fisher_df$`Case wo QV`) * 100)
  fisher_df$`Ctrl % w/ QV` <- sprintf("%.1f", fisher_df$`Ctrl w QV` / (fisher_df$`Ctrl w QV` + fisher_df$`Ctrl wo QV`) * 100)

  if (show_or <- as.logical(arguments$show_or)) {
    or_p_label <- sprintf("%.1f (%s)", as.numeric(as.character(fisher_df$OR)), p_val_text)
  } else {
    or_p_label <- ""
  }

  if (show_per_w_p <- as.logical(arguments$show_per_w_p)) {
    or_p_label <- sprintf("%.1f (%s), %s%% vs. %s%%", as.numeric(as.character(fisher_df$OR)), p_val_text, fisher_df$`Case % w/ QV`, fisher_df$`Ctrl % w/ QV`)
  }

  if (show_ci_w_p <- as.logical(arguments$show_ci_w_p)) {
    or_p_label <- sprintf("%s, [%.2f %.2f]", or_p_label, as.numeric(as.character(fisher_df$CI_low)), as.numeric(as.character(fisher_df$CI_high)))
  }
  ifelse(show_x_axis_title, axis_title_label <- "Odds Ratio", axis_title_label <- "")
  significance_line_position <- 1
  y_breaks_var <- pub_breaks

  # Plot Figure----
  pub_figure_obj <- ggplot(data = fisher_df, aes(x = x_val, y = (as.numeric(as.character(OR))), ymin = CI_low, ymax = CI_high)) +
    theme_cowplot(font_size = as.numeric(arguments$cowplot_font_size)) +
    coord_flip() +
    annotate("text",
      x = fisher_df$x_val + p_line_space,
      y = (as.numeric(as.character(fisher_df$OR))), label = or_p_label,
      hjust = p_val_hjust, size = p_val_textsize, alpha = p_alpha
    )
  if (label_w_10_power) {
    pub_figure_obj <- pub_figure_obj +
      scale_y_continuous(limits = c(y_min_pub, y_max_pub), breaks = y_breaks_var, name = axis_title_label, expand = c(0, 0), trans = "log10", labels = trans_format("log10", math_format(10^.x)))
  } else {
    pub_figure_obj <- pub_figure_obj +
      scale_y_continuous(limits = c(y_min_pub, y_max_pub), breaks = y_breaks_var, name = axis_title_label, expand = c(0, 0), trans = "log10")
  }

  vertical_axis_labels <- fisher_df$Publication_Name

  pub_figure_obj <- pub_figure_obj +
    geom_pointrange(aes(col = eval(parse(text = color_var))), size = as.numeric(arguments$pointrange_size), shape = 15) + # Group
    geom_hline(aes(yintercept = significance_line_position), size = 0.5, alpha = as.numeric(arguments$sig_line_alpha)) + # significance line veritcal
    scale_colour_manual(values = palette_colors) +
    scale_x_continuous(limits = c(0 + y_bottom_buffer, max(fisher_df$x_val) + y_top_buffer), expand = c(0, 0), labels = vertical_axis_labels, breaks = fisher_df$x_val) +
    theme(
      plot.margin = unit(margin_vector, "cm"), axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank()
    )
  # Add legend. Accepts positions outside and inside plot
  if (arguments$pub_leg_pos %in% c("none", "right", "bottom", "top")) {
    pub_figure_obj <- pub_figure_obj + theme(legend.position = arguments$pub_leg_pos)
  } else {
    leg_pos_vector <- as.numeric(unlist(strsplit(arguments$pub_leg_pos, ",")))
    pub_figure_obj <- pub_figure_obj + theme(legend.position = c(leg_pos_vector[1], leg_pos_vector[2]))
  }

  ggsave(plot = pub_figure_obj, filename = paste(output_directory, time_case_prefix, "publication.pdf", sep = ""), width = pub_width, height = pub_height, units = "in", dpi = 300)

  logr::log_print("Finished forest plot")
  logr::log_print("starting accompanying table 1")

  # accompanying table----

  table_for_plot <- fisher_df %>% select(Model = Publication_Name, `Case % w/ QV`, `Ctrl % w/ QV`)
  table_for_plot$`OR (95% CI)` <- sprintf("%.2f (%.2f-%.2f)", as.numeric(as.character(fisher_df$OR)), as.numeric(as.character(fisher_df$CI_low)), as.numeric(as.character(fisher_df$CI_high)))
  logr::log_print(table_for_plot)

  table_figure_left <- pub_figure_obj
  table_figure_left$layers <- NULL
  start_vector <- c(0.6, 4.25, 6.1, 7.65)
  table_figure_left <- table_figure_left +
    theme_cowplot(font_size = 12) + # , font_family = font_variable + scale_color_manual(values=palette_colors) +
    annotate("text", x = fisher_df$x_val, y = rep(start_vector[1], nrow(fisher_df)), label = table_for_plot$Model, size = p_val_textsize, hjust = 0) +
    annotate("text", x = fisher_df$x_val, y = rep(start_vector[2], nrow(fisher_df)), label = table_for_plot$`OR (95% CI)`, size = p_val_textsize, hjust = 0) +
    annotate("text", x = fisher_df$x_val, y = rep(start_vector[3], nrow(fisher_df)), label = table_for_plot$`Case % w/ QV`, size = p_val_textsize, hjust = 0) +
    annotate("text", x = fisher_df$x_val, y = rep(start_vector[4], nrow(fisher_df)), label = table_for_plot$`Ctrl % w/ QV`, size = p_val_textsize, hjust = 0) +
    annotate("text", x = max(fisher_df$x_val) + 0.25, y = start_vector, label = names(table_for_plot)[c(1, 4, 2, 3)], size = p_val_textsize, hjust = 0) +
    scale_y_continuous(limits = c(.5, 9), expand = c(0, 0)) +
    theme(
      axis.line.y = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line.x = element_line(color = "white"), axis.title.x = element_text(color = "white"), axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white")
    )
  ggsave(plot = table_figure_left, filename = paste(output_directory, time_case_prefix, "publication_table_1.pdf", sep = ""), width = pub_width, height = pub_height, units = "in", dpi = 300)

  table_figure_right <- pub_figure_obj
  table_figure_right$layers <- NULL
  start_vector <- c(0.6, 4.25, 6.1, 7.65)
  table_figure_right <- table_figure_right +
    theme_cowplot(font_size = 12) + # , font_family = font_variable + scale_color_manual(values=palette_colors) +
    annotate("text", x = fisher_df$x_val, y = rep(0.2, nrow(fisher_df)), label = p_val_text, size = p_val_textsize, hjust = 0) +
    annotate("text", x = max(fisher_df$x_val) + 0.25, y = 0.2, label = "P", size = p_val_textsize, hjust = 0, fontface = "italic") +
    scale_y_continuous(limits = c(0.1, 1.1), expand = c(0, 0)) +
    theme(
      axis.line.y = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line.x = element_line(color = "white"), axis.title.x = element_text(color = "white"), axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white")
    )
  ggsave(plot = table_figure_right, filename = paste(output_directory, time_case_prefix, "publication_table_2.pdf", sep = ""), width = pub_width, height = pub_height, units = "in", dpi = 300)



  temp <- cowplot::plot_grid(plotlist = list(table_figure_left, pub_figure_obj, table_figure_right), nrow = 1, rel_widths = rel_width_vector) # , rel_widths = c(1,0.5)
  temp <- temp + theme(text = element_text(family = "Arial"))
  ggsave(plot = temp, paste0(output_directory, time_case_prefix, "table_figure.pdf", sep = ""), width = pub_width * 2, height = pub_height, units = "in", dpi = 300)
  # ggsave(plot = temp, paste0(output_directory,time_case_prefix,"table_figure.eps",sep = ""),width = pub_width*2, height = pub_height, units = "in", dpi = 300)
}

#' Initialize log file
#' @param arguments The arguments list with the scripts arguments
initialize_logfile <- function(arguments) {
  dir.create(here("Data"))
  dir.create(log_folder <- here("Data", "fp_forest_plotLog"))
  time_case_prefix <- paste0(gsub(":", "_", gsub(" ", "_", gsub("-", "_", Sys.time()))), "_", arguments$title_str, "_")
  logr::log_open(paste0(log_folder, "/", time_case_prefix, "fp_forest_plot_logfile.log"))
  logr::log_print(arguments)
}


#' Generate statistics from model
#' @param arguments The arguments list with the scripts arguments
#' @param fp_models The arguments list with the scripts arguments
generate_stats <- function(arguments, fp_models) {
  logr::log_print("Entering generate_stats")
  # fp_models[1,] <- fp_models

  model_data <- mclapply(1:nrow(fp_models), function(x) cluster_sample_and_genotype_list_fp(fp_models$path_to_git[x], fp_models$resolution[x], fp_models$min_sample[x], fp_models$Model[x], matrix_geno = "matrix"),
    mc.cores = as.numeric(arguments$cores)
  )
  stats_df <- mclapply(model_data, function(x) atav_fisher_fp(x$genotype_df, x$sample_df, c("case", "ctrl"), cmh_test = TRUE),
    mc.cores = as.numeric(arguments$cores)
  )
  names(stats_df) <- fp_models$Analysis
  logr::log_print("Leaving generate_stats")
  return(stats_df)
}

#' converts the statistical output into a df for use for creating a forest plot
#' @param stats_df statistical dataframe
#' @param fp_models a data frame holding the csv that has the model information
convert_stats_to_fp_obj <- function(stats_df, fp_models) {
  logr::log_print("Entering convert_stats_to_fp_obj")
  OR <- sapply(stats_df, function(x) x[[1]]$estimate)
  CI_low <- sapply(stats_df, function(x) x[[1]]$conf.int[1])
  CI_high <- sapply(stats_df, function(x) x[[1]]$conf.int[2])
  P_uncorrected <- sapply(stats_df, function(x) x[[1]]$p.value)
  case_w_qv <- sapply(stats_df, function(x) x[[2]][1])
  case_wo_qv <- sapply(stats_df, function(x) x[[2]][3])
  ctrl_w_qv <- sapply(stats_df, function(x) x[[2]][2])
  ctrl_wo_qv <- sapply(stats_df, function(x) x[[2]][4])

  logr::log_print("Leaving convert_stats_to_fp_obj")
  return(fp_obj <- data.frame(
    OR = OR, CI_low = CI_low, CI_high = CI_high,
    Analysis = names(stats_df), "Case w QV" = case_w_qv,
    "Case wo QV" = case_wo_qv,
    "Ctrl w QV" = ctrl_w_qv,
    "Ctrl wo QV" = ctrl_wo_qv, P_uncorrected = P_uncorrected,
    P_corrected = P_uncorrected, Publication_Name = fp_models$row_name,
    check.names = FALSE
  ))
}

tryCatch(initialize_logfile(arguments), error = function(e) "Error initializing logfile")
logr::log_print("Main")
if (arguments$path_to_fp_models == "NA") {
  logr::log_print("arguments$path_to_fp_models must be NA")
  tryCatch(display_fp(arguments), error = function(e) stop("Error displaying fp"))
} else {
  fp_models <- fread(arguments$path_to_fp_models) %>% map_df(rev)
  withCallingHandlers(stats_df <- generate_stats(arguments, fp_models), error = function(e) {
    message("error:\n", e)
  })
  fp_obj <- convert_stats_to_fp_obj(stats_df, fp_models)
  # tryCatch(fp_obj <- convert_stats_to_fp_obj(stats_df), error = function (e)  "Error generating fp_obj")
  arguments$fp_obj <- fp_obj
  display_fp(arguments)
  # tryCatch(display_fp(arguments), error = function (e)  "Error displaying fp")
}

logr::log_close()
