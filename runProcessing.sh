#!/bin/bash

Rscript src/processing/0.cleanMethylation-script.R
Rscript src/processing/1.cleanExpression-script.R
./src/processing/runPooled.sh
