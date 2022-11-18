if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager",repos="https://cloud.r-project.org/")

BiocManager::install("biomaRt")
BiocManager::install("GenomicDataCommons")
BiocManager::install("TCGAbiolinks")

# install packages from CRAN
install.packages(c("data.table",
                   "dplyr",
                   "ggplot2",
                   "ggpubr",
                   "gridExtra",
                   "Hmisc",
                   "huge",
                   "igraph",
                   "magrittr",
                   "pROC",
                   "stringr",
                   "tidyverse",
                   "xtable"),
		   repos="https://cloud.r-project.org/")


BiocManager::install("fgsea")