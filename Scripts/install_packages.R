# install required packages
pkgs = c("data.table", "tidyverse", "broom", "openxlsx", "readxl", "here", "umap", "Seurat", "Matrix", "ggrepel", "gridExtra", "gtable", "grid", "parallel", "optparse","docopt","logr","cowplot","ggpubr","stringr","sjPlot","yaml","qqman","GenABEL","broom")
pkgs.na = pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(pkgs.na) > 0) {
  install.packages(pkgs.na, dependencies = TRUE)
}
