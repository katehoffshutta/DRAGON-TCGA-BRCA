#!/bin/bash

# 20220511 Test DRAGON on the four different types: no transformation, npn only, m only, m followed by npn

echo "[DRAGON] processing untransformed data"
python runDragon.py ../../data/interim/gene_expression_processed.tsv ../../data/interim/methylation_processed_no_transf.tsv untransformed

echo "[DRAGON] processing npn data"
python runDragon.py ../../data/interim/gene_expression_processed.tsv ../../data/interim/methylation_processed.tsv npn

echo "[DRAGON] processing m-value data"
python runDragon.py ../../data/interim/gene_expression_processed.tsv ../../data/interim/methylation_processed_m.tsv mvalue

echo "[DRAGON] processing m-value data"
python runDragon.py ../../data/interim/gene_expression_processed.tsv ../../data/interim/methylation_processed_m_npn.tsv mvalue_npn
