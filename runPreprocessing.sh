#!/bin/bash

Rscript src/preprocessing/0.downloadMethylation.R
Rscript src/preprocessing/1.downloadGeneExpression.R
Rscript src/preprocessing/2.mapProbes.R
Rscript src/preprocessing/3.calculateMeanMethylation-script.R
Rscript src/preprocessing/4.calculateGeneExpression.R
Rscript src/preprocessing/5.getClinicalData.R
Rscript src/preprocessing/6.getSubtypes.R
Rscript src/preprocessing/7.makeAnalysisDataset.R
