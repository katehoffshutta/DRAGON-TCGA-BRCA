library(data.table)
library(GenomicDataCommons)


ge_manifest = files() %>%
    GenomicDataCommons::filter( cases.project.project_id == 'TCGA-BRCA') %>% 
    GenomicDataCommons::filter( type == 'gene_expression' ) %>%
    GenomicDataCommons::filter( data_type == 'Gene Expression Quantification') %>%
    GenomicDataCommons::filter( cases.demographic.gender == 'female') %>%
    manifest() %>%
    dplyr::filter(!grepl("splice_junctions",file_name))

tfList = read.table("data/external/TF_names_v_1.01.txt")[,1]

# There are 1639 genes in the TF list. Read a first file to populate the data frame.

direname = ge_manifest$id[1]
filename = ge_manifest$file_name[1]

exp_raw = fread(paste(c("data/external/tcga_BRCA_gene_expression",direname,filename),collapse="/"),sep="\t",header=T)

# remove a few metadata
expression = exp_raw[-c(1:4),] 
tpms_df = expression %>% dplyr::filter(gene_name %in% tfList) %>%
  dplyr::select(gene_id,gene_name)
raw_df = tpms_df
rm(expression)


# Iterate through all of the files
if(!dir.exists("data/interim/")) dir.create("data/interim/")
for(j in 1:length(ge_manifest$id))
{
  print(paste("Processing expression file: ",j))
  # get expression file from s3
  direname = ge_manifest$id[j]
  filename = ge_manifest$file_name[j]

  exp_raw = fread(paste(c("data/external/tcga_BRCA_gene_expression",direname,filename),collapse="/"),sep="\t",header=T)
  

  # remove a few metadata
  expression = exp_raw[-c(1:4),] 

  # filter for TFs
  this_expression = expression %>% dplyr::filter(gene_name %in% tfList) %>% dplyr::select(gene_id,tpm_unstranded)
  
  # merge on gene_id to allow for splice variants
  tpms_df = merge(tpms_df,this_expression,by="gene_id",all.x=T,all.y=T)
  
  # add the sample label
  names(tpms_df)[j+2]=ge_manifest$id[j]
  write.table(tpms_df,file="data/interim/gene_expression.tsv",sep="\t",row.names=F,quote=F)
}


# Perform again selecting the raw counts so that coverage can be assessed in the case of duplicate samples:

for(j in 1:length(ge_manifest$id))
{
  print(paste("Processing expression file: ",j))
  # get expression file from s3
  direname = ge_manifest$id[j]
  filename = ge_manifest$file_name[j]
  
  exp_raw = fread(paste(c("data/external/tcga_BRCA_gene_expression",direname,filename),collapse="/"),sep="\t",header=T)

  # remove a few metadata
  expression = exp_raw[-c(1:4),] 

  # filter for TFs
  this_expression = expression %>% dplyr::filter(gene_name %in% tfList) %>% dplyr::select(gene_id,unstranded)
  
  # merge on gene_id to allow for splice variants
  raw_df = merge(raw_df,this_expression,by="gene_id",all.x=T,all.y=T)
  
  # add the sample label
  names(raw_df)[j+2]=ge_manifest$id[j]
  write.table(raw_df,file="data/interim/gene_expression_unstranded_raw.tsv",sep="\t",row.names=F,quote=F)
}
