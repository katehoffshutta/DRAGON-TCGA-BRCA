library(data.table)
library(tidyverse)

# The purpose of this code is map probes to gene, going from a matrix of probes x 
# samples to a matrix of genes x samples. 
# At this state, we want to simply map each promoter to its nearest gene. 
# For this we need the probe manifest. 
# See this helpful file: https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/infinium-methylationepic-manifest-column-headings.pdf -->


# source hg38 with gencode 36 from https://zwdzwd.github.io/InfiniumAnnotation
download.file('https://zhouserver.research.chop.edu/InfiniumAnnotation/20210615/HM450/HM450.hg38.manifest.gencode.v36.tsv.gz', destfile = "data/external/HM450.hg38.manifest.gencode.v36.tsv.gz")

# unzip
system2(command="gunzip",args=c("data/external/HM450.hg38.manifest.gencode.v36.tsv.gz"))

# load into memory
manifest = data.frame(fread("data/external/HM450.hg38.manifest.gencode.v36.tsv",sep="\t",header=T))

# Order of operations: for each probe, need to (1) see if there is >1 gene within TSS-1500 (2) if so, see if they all correspond to the same gene (splice isoforms?), (3) if so, label probe with gene, (4) if not, map to the closer gene. Note this is arbitrary choice!! May make sense to allow to map to both genes?

processRow = function(x) # x is one row of the manifest
{
  x = as.data.frame(x)
  genes = str_split(x$geneNames,";",simplify=T)
  tssDist = as.numeric(str_split(x$distToTSS,";",simplify=T))
  promoterRegions = which(tssDist > -1500 & tssDist < 0) # negative = upstream of TSS
  promoterMap = genes[promoterRegions]
  if(length(unique(promoterMap)) > 1)
  {
    index = which.max(tssDist[promoterRegions])
    thisRow = c(x$probeID,promoterMap[index],tssDist[promoterRegions][index])
    return(thisRow)
  }
  
  thisGene = ifelse(length(promoterMap)>0,unique(promoterMap), NA)
  thisDist = ifelse(length(promoterMap)>0,max(tssDist[promoterRegions]),NA) 
  return(c(x$probeID,thisGene,thisDist))
}

# perform the mapping

mymap = matrix(rep(NA,3*nrow(manifest)),ncol=3)
for(i in 1:nrow(manifest))
{
  if(i %% 10000 == 0) print(paste("probe number:",i))
  mymap[i,]=processRow(manifest[i,])
}

# sanity check that nothing is mapped to multiple genes

grep(";",mymap[,2])

# Write table for those probes that mapped to promoter regions:
if(!dir.exists("data/processed")) dir.create("data/processed")

write.table(na.omit(mymap),file="data/processed/PromoterProbeMap450K.txt",sep="\t",row.names=F,quote=F)
