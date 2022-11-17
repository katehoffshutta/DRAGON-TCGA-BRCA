#!/bin/bash

Rscript src/postprocessing/0.full-postprocess.R
Rscript src/postprocessing/1.tf_nbhds.R
Rscript src/postprocessing/2.ZNF332_subtype_boxplots.R
Rscript src/postprocessing/3a.community5_analysis.R
Rscript src/postprocessing/3b.TFAP2B_classifier.R
Rscript src/postprocessing/4.community38_analysis.R

# left out ma analysis for now
