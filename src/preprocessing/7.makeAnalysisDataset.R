# The purpose of this script is to combine methylation data, UUIDs, TCGA short and long IDs, batch variables such as sample, plate, vial, etc., clinical and demoraphic data, into one dataset that will be frozen for analysis.

# columns: 
# UUID
# TCGA long
# TCGA short
# sample
# vial
# portion
# analyte
# plate
# center
# sex
# mrna subtype
# sample subtype 
# gene1_meth_mean
# ...
# geneN_meth_mean
# gene1_meth_sd
# ...
# geneN_meth_sd

library(data.table)
library(dplyr)
library(stringr)
library(TCGAbiolinks)

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

subtypeDat = read.table("data/external/brcaTypes.tsv",sep="\t",header=T)
analysisDataset_B = data.frame("TCGA_long"=subtypeDat$pan.samplesID, 
                               "TCGA_short"=subtypeDat$shortID,
                               "TCGA_shorter"=substr(subtypeDat$shortID,start=1,stop=12),
                               "Subtype_mRNA"=subtypeDat$Subtype_mRNA)

dim(analysisDataset_B)
# [1] 1218    4

demoDat = read.table("data/external/TCGA_BRCA_clinical.tsv",sep = "\t", header=T)
analysisDataset_C = data.frame("TCGA_shorter" = substr(demoDat$submitter_id,start=1,stop=12),
                               "race" = demoDat$race,
                               "ethnicity" = demoDat$ethnicity,
                               "gender" = demoDat$gender)

dim(analysisDataset_C)
# [1] 1085    4

methylationAll = read.table("data/interim/beta_means.txt")
dim(methylationAll)
#[1]  885 1591

names(methylationAll) = paste0(names(methylationAll),"_methylation")
#names(methylationAll)[1]="UUID"

#UUIDs = methylationAll[,1]%>% 
#  str_replace("X","") %>%
#  str_replace_all("\\.","-")

analysisDataset_D = methylationAll 
analysisDataset_D$UUID = row.names(methylationAll)
dim(analysisDataset_D)
#[1]  885 1591
# 885 is the right number for the analysis dataset

results = GDCquery(project = "TCGA-BRCA",
                   data.category = "DNA Methylation",
                   platform = "Illumina Human Methylation 450",
                   legacy=FALSE)

resultsDF = getResults(results) %>% dplyr::filter(data_format == "TXT")

analysisDataset_E = data.frame("UUID" = resultsDF$id,
                               "TCGA_short" = resultsDF$sample.submitter_id,
                               "TCGA_long_me_th" = resultsDF$cases,
                               "GDCquery_subtype" = resultsDF$sample_type)

for(i in 1:nrow(analysisDataset_E))
{
  split_tcga = str_split(analysisDataset_E$TCGA_long_me_th[i],"-",simplify=T)
  analysisDataset_E$sample[i] = substr(split_tcga[1,4],1,2)
  analysisDataset_E$vial[i] = substr(split_tcga[1,4],3,3)
  analysisDataset_E$portion[i] = substr(split_tcga[1,5],1,2)
  analysisDataset_E$analyte[i] = substr(split_tcga[1,5],3,3)
  analysisDataset_E$plate[i] = split_tcga[6]
  analysisDataset_E$center[i] = split_tcga[7]
}

# remove anything that is normal tissue
analysisDataset_E_tumor = analysisDataset_E %>% dplyr::filter(sample < 10)

# analysisDataset_BE is identifiers and subtype
analysisDataset_EB = left_join(analysisDataset_E_tumor, analysisDataset_B, by = "TCGA_short")
dim(analysisDataset_EB)
# [1] 798  13

# join the clinical data to get self-reported gender 
# intentionally noting this is not necessarily sex
# race/ethnicity
analysisDataset_BEC = left_join(analysisDataset_EB, analysisDataset_C, by = "TCGA_shorter") %>% 
  dplyr::filter(gender == "female")
apply(analysisDataset_BEC,2,function(x){summary(as.factor(x))})
dim(analysisDataset_BEC)
# [1] 778 16

# join the  methylation data with the GDC query to map UUIDs to TCGA_short
analysisDataset_DE = left_join(analysisDataset_D, analysisDataset_E,by = "UUID") %>% 
  dplyr::select(-c(TCGA_long_me_th,GDCquery_subtype,sample,vial,portion,analyte,plate,center))
dim(analysisDataset_DE)
#[1]  885 1592

problems = analysisDataset_DE %>% dplyr::filter(TCGA_short %in% c("TCGA-A7-A26E-01A","TCGA-A7-A26J-01A"))
methIndx = grep("_methylation",names(problems))
cor(as.numeric(problems[1,methIndx]),as.numeric(problems[4,methIndx]),use="complete.obs")
cor(as.numeric(problems[2,methIndx]),as.numeric(problems[3,methIndx]),use="complete.obs")

# the correlation is high, so omit two at random. 
set.seed(42)
index1 = sample(c(1,4),size=1)
index2 = sample(c(2,3),size=1)
problemUUIDs = problems[c(index1,index2),]$UUID
analysisDataset_DE_uniq = analysisDataset_DE %>% dplyr::filter(!(UUID %in% problemUUIDs))
dim(analysisDataset_DE_uniq)
# [1]  883 1592
analysisDataset_BEC_uniq = analysisDataset_BEC %>% dplyr::filter(!(UUID %in% problemUUIDs))
dim(analysisDataset_BEC_uniq)
# [1] 776  16

expr = read.table("data/interim/gene_expression.tsv",sep="\t",header=T)
names(expr)[1]="ens_id"
names(expr)[2]="gene_name"

# remove any genes with the same value across all data
means = apply(expr[,-c(1:2)],1,mean)
sds = apply(expr[,-c(1:2)],1,sd)
means[which(sds==0)]

# we have only one gene that is mapped to two ensembl IDs and it's a little strange
# ENSG00000228623.6    ZNF883 
# ENSG00000285447.1    ZNF883

# Exclude this TF for now and ask about it later
exprClean = expr[-union(which(sds==0),which(expr$gene_name=="ZNF883")),]
length(unique(exprClean$gene_name))
length(exprClean$gene_name)

# now we can wrangle in to the same form as methylation: column = gene_name, row = UUID
gene_names = exprClean$gene_name
sample_names = colnames(exprClean)[-c(1:2)]
exprT = data.frame(t(exprClean[,-c(1:2)]))
colnames(exprT) = paste0(gene_names,"_expr")

row.names(exprT) = row.names(exprT)%>% 
  str_replace("X","") %>%
  str_replace_all("\\.","-")
  
  
exprTCGA = TCGAtranslateID(row.names(exprT))
names(exprTCGA)[1]="UUID"
names(exprTCGA)[2]="TCGA_short"

exprT$UUID = rownames(exprT)
analysisDataset_F = merge(exprTCGA,exprT, by="UUID")
# problems are TCGA-A7-A0DB-01A and TCGA-A7-A13D-01A and TCGA-A7-A13E-01A and
# TCGA-A7-A26E-01A and TCGA-A7-A26J-01A
problems1 = analysisDataset_F %>% dplyr::filter(TCGA_short == c("TCGA-A7-A13D-01A"))
problems2 = analysisDataset_F %>% dplyr::filter(TCGA_short == c("TCGA-A7-A13E-01A"))
problems3 = analysisDataset_F %>% dplyr::filter(TCGA_short == c("TCGA-A7-A26E-01A"))
problems4 = analysisDataset_F %>% dplyr::filter(TCGA_short == c("TCGA-A7-A26J-01A"))

expr_raw = read.table("data/interim/gene_expression_unstranded_raw.tsv",sep="\t",header=T)
colnames(expr_raw) = colnames(expr_raw)%>% 
  str_replace("X","") %>%
  str_replace_all("\\.","-")
  
# calculate read depth by summing counts across all genes; choose the expt w best depth
problemUUIDs = analysisDataset_F %>% dplyr::filter(TCGA_short %in% c("TCGA-A7-A13D-01A",
                                         "TCGA-A7-A13E-01A",
                                         "TCGA-A7-A26E-01A",
                                         "TCGA-A7-A26J-01A")) %>%
  dplyr::select(UUID)

tfList = read.table("data/external/TF_names_v_1.01.txt")[,1]

read_depth = expr_raw %>% dplyr::filter(gene_name %in% tfList) %>%
  dplyr::select(all_of(problemUUIDs$UUID)) %>% 
  apply(2,sum) %>% data.frame()

read_df = data.frame("UUID"=row.names(read_depth),"read_depth"=read_depth[,1])
problems1_depth = merge(problems1,read_df,by=c("UUID"))
problems2_depth = merge(problems2,read_df,by=c("UUID"))
problems3_depth = merge(problems3,read_df,by=c("UUID"))
problems4_depth = merge(problems4,read_df,by=c("UUID"))

exclude1 = problems1[which.min(problems1_depth$read_depth),]$UUID
exclude2 = problems2[which.min(problems2_depth$read_depth),]$UUID
exclude3 = problems3[which.min(problems3_depth$read_depth),]$UUID
exclude4 = problems4[which.min(problems4_depth$read_depth),]$UUID

analysisDataset_F_uniq = analysisDataset_F %>% dplyr::filter(!UUID %in% c(exclude1,exclude2,exclude3,exclude4))

# left join the GDC query/subtype/clinical data with the methylation/GDC query data
analysisDataset_BECD = inner_join(analysisDataset_BEC_uniq,analysisDataset_DE_uniq,by=c("UUID","TCGA_short"))
dim(analysisDataset_BECD)
# [1]  776 1606

# inner join: take stuff that has both methylation and expression
analysisDataset_BECDF = inner_join(analysisDataset_BECD, analysisDataset_F_uniq,by=c("TCGA_short")) %>%
  dplyr::rename(UUID_meth = UUID.x, UUID_exp = UUID.y)
dim(analysisDataset_BECDF)
# [1]  770 3236

# Is there anything weird left over?
table1data = analysisDataset_BECDF %>% dplyr::select(-contains("_expr")) %>% dplyr::select(-contains("_methylation"))
apply(table1data,2,function(x){summary(as.factor(x))})

# Note that there are four samples that have both a primary tumor and a metastasis
table1data %>% dplyr::filter(TCGA_shorter %in% c("TCGA-AC-A6IX",
                                                 "TCGA-BH-A1ES",
                                                 "TCGA-BH-A1FE",
                                                 "TCGA-E2-A15K")) %>%
  arrange(TCGA_shorter)

# remove metastases, since there are only 5
analysisDataset_final = analysisDataset_BECDF %>% dplyr::filter(GDCquery_subtype != "Metastatic")
dim(analysisDataset_final)
#[1]  765 3236

table1final = analysisDataset_final%>% dplyr::select(-contains("_expr")) %>% dplyr::select(-contains("_methylation")) %>%
  dplyr::select(-contains("UUID")) %>% dplyr::select(-contains("TCGA"))

apply(table1final,2,function(x){summary(as.factor(x))})
length(grep("_expr",names(analysisDataset_final)))
length(grep("_meth",names(analysisDataset_final)))

write.table(analysisDataset_final,file="data/interim/analysis_dataset.tsv",sep="\t",row.names=F,quote=F)
