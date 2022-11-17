if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
BiocManager::install("GenomicDataCommons")
BiocManager::install("TCGAbiolinks")

# install packages from CRAN
install.packages(c("data.table",
                   "dplyr",
                   "fgsea",
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
                   "xtable"))
