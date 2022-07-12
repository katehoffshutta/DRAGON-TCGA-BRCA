# Get TCGA ids, merge with clinical data, select only 450k values, impute missing values, log transform, standardize

library(aws.s3)
library(data.table)
library(dplyr)
library(GenomicDataCommons)
library(magrittr)
library(huge)

# set the appropriate profile
Sys.setenv("AWS_PROFILE" = "MFA")

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

meanImpute = function(x)
{
  x[is.na(x)] <- mean(x,na.rm=T)
  return(x)
}

# this function removes genes that weren't measured anywhere and should be done agnostic to subtype
removeMissing = function(meth_df, thres = 0.2)
{
  propMiss = apply(meth_df,2,function(x){sum(is.na(x))/length(x)})
  skipIndices = which(propMiss > thres)
  print(paste("Number of genes with > 0.2 missing:",length(skipIndices)))
  print("Genes omitted:")
  print(names(meth_df)[skipIndices])
  meth_dfClean = meth_df[,-skipIndices]
  return(meth_dfClean)
}

betaToM = function(beta)
{
  return(log2(beta/(1-beta))) # reference: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
}

# this function does imputation and transformation and should be applied within subtype
cleanMethylationData = function(meth_df, npn=T, mval=F) # meth_df is a data frame of beta means, rows=samples, first column=sample ids,  cols=genes
{
  meth_dfComplete = removeMissing(meth_df,thres = 0.2)
  # skip anything w/more than 0.2 missing
  
  for(i in 2:ncol(meth_dfComplete))
  {
    meth_dfComplete[,i] = meanImpute(meth_df[,i])
  }
  
  # sanity check the imputation
  summary(apply(meth_dfComplete,2,function(x){sum(is.na(x))}))
  summary(apply(meth_dfComplete[,-1],2,sd))
  
  if(mval & !npn)
  {
    transf_data = data.frame(apply(meth_dfComplete[,-1],2,betaToM))
    row.names(transf_data)= meth_dfComplete[,1]
    return(transf_data)
  }
  
  if(npn)
  {
    # do the nonparanormal transformation
    transf_data = data.frame(huge.npn(meth_dfComplete[,-1]))
    row.names(transf_data)= meth_dfComplete[,1]
    return(transf_data)
  }
  
  if(mval & npn)
  {
    # first apply M value transformation
    transf_data_m = apply(meth_dfComplete[,-1],2,betaToM)
    # next apply npn
    transf_data_m_npn = data.frame(huge.npn(transf_data_m))
    row.names(transf_data_m_npn)= meth_dfComplete[,1]
    return(transf_data_m_npn)
  }
  
  row.names(meth_dfComplete) = meth_dfComplete[,1]
  
  return(meth_dfComplete[,-1])
  
}

save_object(object="analysis_dataset.tsv",
                             bucket = "netzoo/supData/dragon/dragonInputData", 
                             region="us-east-2",
                             multipart=F,
                             file="analysis_dataset.tsv")

analysis_dataset = read.table("analysis_dataset.tsv",sep="\t",header=T)

dim(analysis_dataset) 
# [1] 765 3236

transformedData = analysis_dataset %>%
  dplyr::select(TCGA_short,contains("_methylation")) %>%
  cleanMethylationData(npn=T,mval=F)

standardize = function(x){return((x-mean(x))/sd(x))}
methTransfStd = data.frame(apply(transformedData,2,standardize))

summary(apply(methTransfStd,2,mean))
summary(apply(methTransfStd,2,sd))

methTransfStd$TCGA_barcode = row.names(methTransfStd)
outdf = methTransfStd %>% dplyr::relocate(TCGA_barcode)

write.table(outdf,file="../../data/interim/methylation_processed_npn.tsv",sep="\t",row.names=F,quote=F)
put_object(file="../../data/interim/methylation_processed_npn.tsv",
           bucket = "netzoo/supData/dragon/dragonInputData", 
           region="us-east-2",
           multipart=F)

# for sensitivity analyses, also store the untransformed data

rawData = analysis_dataset %>%
  dplyr::select(TCGA_short,contains("_methylation")) %>%
  cleanMethylationData(npn=F,mval=F)
rawData$TCGA_barcode = row.names(rawData)

outdf = rawData %>% dplyr::relocate(TCGA_barcode)
write.table(outdf,file="../../data/interim/methylation_raw.tsv",sep="\t",row.names=F,quote=F)
put_object(file="../../data/interim/methylation_raw.tsv",
           bucket = "netzoo/supData/dragon/dragonInputData", 
           region="us-east-2",
           multipart=F)

# for sensitivity analyses, also store the M-value transformation

mvalData = analysis_dataset %>%
  dplyr::select(TCGA_short,contains("_methylation")) %>%
  cleanMethylationData(npn=F,mval=T)
mvalData$TCGA_barcode = row.names(rawData)
outdf = mvalData %>% dplyr::relocate(TCGA_barcode)

write.table(outdf,file="../../data/interim/methylation_processed_m_val.tsv",sep="\t",row.names=F,quote=F)
put_object(file="../../data/interim/methylation_processed_m_val.tsv",
           bucket = "netzoo/supData/dragon/dragonInputData", 
           region="us-east-2",
           multipart=F)
