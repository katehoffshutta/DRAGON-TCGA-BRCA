# get molecular subtypes from TCGA via TCGAbiolinks
# the source of this data is https://www.synapse.org/#!Synapse:syn8402849
# read more at https://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/subtypes.html

library(dplyr)
library(TCGAbiolinks)
library(tidyverse)

molecular.subtypes <- PanCancerAtlas_subtypes()
brcaTypes = molecular.subtypes %>% dplyr::filter(cancer.type == "BRCA") 
brcaTypes$shortID = substr(brcaTypes$pan.samplesID,start=1,stop=16)
write.table(brcaTypes,"../../data/external/brcaTypes.tsv",sep="\t",row.names=F,quote=F)


