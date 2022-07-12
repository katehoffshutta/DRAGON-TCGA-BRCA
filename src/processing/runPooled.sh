#!/bin/bash

echo "[DRAGON] processing pooled data"
python 2.runDragon.py ../../data/interim/gene_expression_processed.tsv ../../data/interim/methylation_processed_npn.tsv ../../data/processed/Pooled Pooled

