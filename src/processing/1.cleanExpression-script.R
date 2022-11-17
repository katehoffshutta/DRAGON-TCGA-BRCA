library(data.table)
library(dplyr)
library(GenomicDataCommons)
library(magrittr)
# function from: https://seandavi.github.io/post/2017-12-29-genomicdatacommons-id-mapping/

TCGAtranslateID = function(file_ids, legacy = FALSE) {
  info = files(legacy = legacy) %>%
    GenomicDataCommons::filter( ~ file_id %in% file_ids) %>%
    GenomicDataCommons::select('cases.samples.submitter_id') %>%
    results_all()
  # The mess of code below is to extract TCGA barcodes
  # id_list will contain a list (one item for each file_id)
  # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
  id_list = lapply(info$cases,function(a) {
    a[[1]][[1]][[1]]})
  # so we can later expand to a data.frame of the right size
  barcodes_per_file = sapply(id_list,length)
  # And build the data.frame
  return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                    submitter_id = unlist(id_list)))
}

analysis_dataset = read.table("data/interim/analysis_dataset.tsv",sep="\t",header=T)

expr = analysis_dataset %>% dplyr::select(TCGA_short,contains("_expr")) 
row.names(expr) = analysis_dataset$TCGA_short

# Handle NAs

nNA = apply(expr,2,function(x){sum(is.na(x))})
summary(nNA) # sanity check, should be zero

# write function for gene count across samples
# omit anything that is not expressed w count >1 in at least 20% of samples
omit = ifelse(apply(expr[,-1],2,function(x){(sum(x>1)/nrow(expr))>0.2}),F,T)

omittedTFs = expr %>% dplyr::select(-TCGA_short) %>%
  dplyr::select(which(omit==T))

exprClean = expr %>% dplyr::select(-TCGA_short) %>%
  dplyr::select(which(omit==F))
# sanity checks
dim(exprClean)
# should all be true
summary(apply(exprClean,2,function(x){sum(x)/nrow(exprClean)>0.2}))

# These are tpm so I think they have already been normalized but should look this up on TCGA. 
# Log transform with pseudocounts: f(x) = log(x+1)

ltpc = function(x){return(log(x+1))}
exprTransf = data.frame(apply(exprClean,2,ltpc))

# standardize 
standardize = function(x){return((x-mean(x))/sd(x))}
exprTransfStd = data.frame(apply(exprTransf,2,standardize))

# sanity check standardization
summary(apply(exprTransfStd,2,mean))
summary(apply(exprTransfStd,2,sd))

exprTransfStd$TCGA_barcode = row.names(exprTransfStd)
outdf = exprTransfStd %>% dplyr::relocate(TCGA_barcode)

# write analysis dataset to file
write.table(outdf,file="data/interim/gene_expression_processed.tsv",sep="\t",row.names=F,quote=F)


