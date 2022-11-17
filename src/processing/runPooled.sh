#!/bin/bash

echo "[DRAGON] processing pooled data"
python src/processing/2.runDragon.py data/interim/gene_expression_processed.tsv data/interim/methylation_processed_npn.tsv Pooled5

