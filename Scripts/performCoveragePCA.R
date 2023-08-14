# performCoveragePCA.R: perfomrs Coverage PCA on a per cluster level

'Usage: 
  performCoveragePCA.R --ready_for_pca_file=<ready_for_pca_file> [--cluster=<cluster>]
  
  Options:
  -h --help
  --ready_for_pca_file=<ready_for_pca_file> The csv.READYforPCA.txt input file
  --cluster=<cluster> current cluster running

' -> doc

library(docopt)
library(devtools)
library(here)
# biplot(pca)0
# library(ggbiplot)
source('https://raw.githubusercontent.com/vqv/ggbiplot/master/R/ggbiplot.r')
source('https://raw.githubusercontent.com/vqv/ggbiplot/master/R/ggscreeplot.r')
arguments <- docopt(doc, version = 'performCoveragePCA.R 1.0')

# debugging
# arguments <- list()
# arguments$sample_file <- ""
# arguments$ready_for_pca_file <- ""

log_folder <- "Data/covPCA_log"
time_case_prefix <- paste0(gsub(":","_",  gsub(" ","_", gsub("-","_",Sys.time()))), "_")
logr::log_open(paste0(log_folder,"/", time_case_prefix, "_", arguments$cluster ,"_covPCA_logfile.log"))

logr::log_print(arguments)
pc_path <- here(paste0("Results/KinshipFlashPCA/"))
logr::log_print(sampleFile <- list.files(pc_path, "*_kinship_pruned_sample.txt$", full.names = TRUE)[1])

file <- arguments$ready_for_pca_file

# Read files
sampleCov <- read.table(file = file, sep = "\t", header = T, stringsAsFactors = F)
samplesDF <- read.table(file = sampleFile, sep = "\t", header = F, stringsAsFactors = F)
# Rename
row.names(sampleCov) <- sampleCov$Samples
sampleT <- sampleCov[, 2:(dim(sampleCov)[2] - 1)] # remove also Type (last column)

# look for constant columns and remove them
sampleT <- sampleT[, apply(sampleT, 2, var, na.rm = TRUE) != 0]

pca <- prcomp(sampleT,
  center = TRUE,
  scale. = TRUE
)

# Visualize the PCAs
for (first in 1:5) {
  pdf(file = paste0(file, "pca", first, "_", (first + 1), ".pdf"))
  g <- ggbiplot(pca,
    choices = c(first, (first + 1)), obs.scale = 1, var.scale = 1,
    groups = as.factor(sampleCov$Sample_TYPE),
    #              labels=rownames(sampleT),
    ellipse = TRUE, var.axes = FALSE,
    circle = TRUE
  )
  g <- g + scale_color_discrete(name = "")
  g <- g + theme(
    legend.direction = "horizontal",
    legend.position = "top"
  )
  print(g)
  dev.off()
}

SD <- 3
toRemove <- ""
# Usually filter up to PC explaining >3%
PCtoFilter <- 6
for (i in 1:PCtoFilter) {
  pc_sd <- sd(pca$x[, i])
  remove <- pca$x[abs(pca$x[, i]) > (pc_sd * SD), ]
  logr::log_print(paste0("Removing ", length(rownames(remove)), " based on PC", i))
  toRemove <- c(toRemove, rownames(remove))
}

samplesToRemove <- unique(toRemove)

afterFilterSampleFile <- samplesDF[!(samplesDF$V1 %in% samplesToRemove), ]
write.table(afterFilterSampleFile,file=paste0(file,"_PCAsd",SD,"_filtered.txt"),quote=FALSE, sep = "\t", qmethod = "double",row.name=F,col.names = F)
