# DRAGON TCGA BRCA Application
This repository presents an illustrative application of the DRAGON algorithm (https://arxiv.org/abs/2104.01690) to promoter methylation and gene expression in breast cancer data from TCGA.

## How to reproduce the analysis
Clone this repository. From the base directory in the repository, run the three shell scripts sequentially:

```
./runPreprocessing.sh
./runProcessing.sh
./runPostprocessing.sh
```

## Workflow
In src/preprocessing, you will find the scripts used to pull TCGA data, to map methylation probes to TF promoter regions, and to merge phenotypic, methylation, and gene expression data.

In src/processing, you will find the scripts used to clean the promoter-level methylation values and the gene expression data along with the scripts to run DRAGON.

In src/postprocessing, you will find the scripts used to analyze the DRAGON results and make the tables and figures for these sections of the paper.
