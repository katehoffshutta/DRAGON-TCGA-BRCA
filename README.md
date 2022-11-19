# DRAGON TCGA BRCA Application

This repository presents an illustrative application of the DRAGON algorithm (https://arxiv.org/abs/2104.01690) to promoter methylation and gene expression in breast cancer data from TCGA. All of the data are downloaded directly from TCGA; running these scripts will require that you have approximately 20G of storage available. 

## How to reproduce the analysis

Clone this repository. To set up the Python and R tools you will need, install and activate the conda environment in `dragon_env.yml` and install R packages using `installPackages.R`. 

From the base directory in the repository, run the three shell scripts sequentially:

```
./runPreprocessing.sh
./runProcessing.sh
./runPostprocessing.sh
```

Reference: [![DOI](https://zenodo.org/badge/513280418.svg)](https://zenodo.org/badge/latestdoi/513280418)

## Workflow
In src/preprocessing, you will find the scripts used to pull TCGA data, to map methylation probes to TF promoter regions, and to merge phenotypic, methylation, and gene expression data.

In src/processing, you will find the scripts used to clean the promoter-level methylation values and the gene expression data along with the scripts to run DRAGON.

In src/postprocessing, you will find the scripts used to analyze the DRAGON results and make the tables and figures for these sections of the paper.

## A note about operating systems

This analysis was run on Linux. On other OS, there is sometimes a difference in GenomicDataCommons functionality: the attribute `file_name` of the GDC manifest may be `filename` (no underscore) and you will need to change this in the scripts.



