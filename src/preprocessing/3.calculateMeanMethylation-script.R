library(data.table)
library(GenomicDataCommons)
library(tidyverse)

probeMap = read.table("data/processed/PromoterProbeMap450K.txt",sep="\t", header = T)
names(probeMap) = c("probeID","gene","distToTSS")
allGenes = unique(probeMap[,2])

# sanity check, should be empty
doubleGenes = allGenes[grep(";",allGenes)] 
doubleGenes

# Calculate average beta values per gene
# from Lambert SA, Jolma A, Campitelli LF, Das PK, Yin Y, Albu M, Chen X, Taipale J, Hughes TR, Weirauch MT.(2018) 
# The Human Transcription Factors. Cell. 172(4):650-665. doi: 10.1016/j.cell.2018.01.029. Review.
# use download here
download.file('http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt', 
              destfile = "data/external/TF_names_v_1.01.txt")

tfList = read.table("data/external/TF_names_v_1.01.txt")[,1]
tfGenes = intersect(allGenes,tfList)

length(tfList)
length(tfGenes)

tfProbes = probeMap[probeMap[,2] %in% tfGenes,1]
length(tfProbes)

# get manifest for only 450k array
ge_manifest = files() %>%
    GenomicDataCommons::filter( cases.project.project_id == 'TCGA-BRCA') %>%
    GenomicDataCommons::filter( type == 'methylation_beta_value' ) %>%
    GenomicDataCommons::filter( platform == "illumina human methylation 450") %>%
    GenomicDataCommons::filter( cases.demographic.gender == 'female') %>%
    manifest()

write.table(ge_manifest,file = "data/external/TCGA_BRCA_methylation_manifest_450k.txt", sep = "\t", row.names = FALSE, quote = FALSE)

head(ge_manifest)

names(probeMap) = c("probeID","gene","distToTSS")

# betas is going to be a n x p matrix where each row is a person and each column is a TF promoter methylation value. 
# we will iterate through people in the outer loop and TFs in the inner loop. 
nSamples = length(ge_manifest$id)
nTFs = length(tfGenes)
betas = matrix(rep(NA,nSamples*nTFs),ncol=nTFs)
betaSDs = matrix(rep(NA,nSamples*nTFs),ncol=nTFs)
betaNs = matrix(rep(NA,1*nTFs),ncol=nTFs)
colnames(betas) = tfGenes
colnames(betaSDs) = tfGenes
colnames(betaNs) = tfGenes

# do one preprocess where we get just the genes
theseProbes = list()
for(i in 1:length(tfGenes))
{
    if(i %% 300 == 0) print(i)
    thisGene = tfGenes[i]
    # we could find the probes separately if we have a bigger list and more efficiency becomes necessary
    theseProbes[[i]] = probeMap %>% dplyr::filter(gene == thisGene) %>% dplyr::select(probeID) %>% unique() %>% dplyr::select(probeID)
}

processFile = function(x)
{
  x = as.data.frame(t(x))
  print(x$id)
  direname = x$id
  filename = x$file_name
  methylation = fread(paste(c("data/external/tcga_BRCA_methylation",direname,filename),collapse="/"))
  # # save_object(object = filename,
  # #             bucket = paste(c("netzoo/supData/dragon/dragonDataPreprocessing/methylation/",direname),collapse=""), 
  # #             region="us-east-2",
  # #             file = "tmp.txt")
  # # 
  # methylation = fread("tmp.txt",sep="\t",header=F)
  names(methylation)=c("probeID","beta")
  # 
  # # remove the file once it's in memory
  # system2(command="rm",args=c("tmp.txt"))
  # 
  # # extract only the probes corresponding to TFs
  tfMethylation = methylation %>% dplyr::filter(probeID %in% tfProbes) 
  return(tfMethylation)
}

a = apply(ge_manifest,1,processFile)
b = a %>% reduce(left_join,by="probeID")
names(b)=c("probeID",ge_manifest$id)
write.csv(b,"data/interim/all_betas.csv")

# extract only the betas we want and calculate their mean
# This will need to be redone after batch correction if there is any to be done.
# No batch correction was needed; see notebooks/pcPlots.Rmd
betaMeans = matrix(NA,nrow = length(tfGenes), ncol = ncol(b)-1)
betaSDs = matrix(NA,nrow = length(tfGenes), ncol = ncol(b)-1)

for(i in 1:length(tfGenes))
{
  if(i %% 300 == 0) print(i)
  thisGene = tfGenes[i]
  theseProbes = probeMap %>% dplyr::filter(gene == thisGene) %>% dplyr::select(probeID) 
  theseBetas = b %>% dplyr::filter(probeID %in% theseProbes$probeID) %>% dplyr::select(-probeID)
  betaMeans[i,] = apply(theseBetas,2,mean,na.rm=T) # any probe that was missed just doesn't contribute to the average
  betaSDs[i,] = apply(theseBetas,2,sd,na.rm=T) 
}

betaMeansDF = data.frame(t(betaMeans))
colnames(betaMeansDF) = tfGenes
row.names(betaMeansDF) = names(b)[-1]
  
betaSDsDF = data.frame(t(betaSDs))
colnames(betaSDsDF) = tfGenes
row.names(betaSDsDF) = names(b)[-1]

# write to local file
write.table(betaMeansDF,file="data/interim/beta_means.txt",quote=T)
write.table(betaSDsDF,file="data/interim/beta_sds.txt",quote=T)


