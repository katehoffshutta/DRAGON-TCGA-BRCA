# libraries
library(aws.s3)
library(data.table)
library(dplyr)
library(fgsea)
library(igraph)
library(tidyverse)

rm(list = ls())
library("biomaRt")
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)

# set AWS profile
Sys.setenv("AWS_PROFILE" = "MFA")

prefix = "Pooled"
thres = 0.005 # arbitrary

# get the results from S3
save_object(object = "dragon_mat.tsv",
            bucket = paste("netzoo/supData/dragon/dragonOutputFiles",prefix,sep="/"),
            region="us-east-2",
            file = paste(c("../../data/processed",prefix,"dragon_mat.tsv"),collapse="/"))

save_object(object = "dragon_adj_p.tsv",
            bucket = paste("netzoo/supData/dragon/dragonOutputFiles",prefix,sep="/"),
            region="us-east-2",
            file = paste(c("../../data/processed",prefix,"dragon_adj_p.tsv"),collapse="/"))

save_object(object = "dragon_input_mat.tsv",
            bucket = paste("netzoo/supData/dragon/dragonOutputFiles",prefix,sep="/"),
            region="us-east-2",
            file = paste(c("../../data/processed",prefix,"dragon_input_mat.tsv"),collapse="/"))


res = data.frame(fread(paste(c("../../data/processed",prefix,"dragon_mat.tsv"),collapse="/"),
                       sep="\t",header=T), row.names = 1)
adj_p_val = data.frame(fread(paste(c("../../data/processed",prefix,"dragon_adj_p.tsv"),collapse="/"),
                         sep="\t",header=T),row.names = 1)
input = data.frame(fread(paste(c("../../data/processed",prefix,"dragon_input_mat.tsv"),collapse="/"),
                         sep="\t",header=T),row.names = 1)

diag(adj_p_val) <- 1

par_cor <- read.table(file = "../OneDrive_1_11/dragon_mat.tsv",sep = "\t", row.names = 1, header = TRUE)
X <- read.table(file = "../OneDrive_1_11/dragon_input_mat.tsv",sep = "\t", row.names = 2, header = TRUE)
X <- X[,-1]
dim(X)
library(Hmisc) # You need to download it first.
cor_X <- rcorr(as.matrix(X), type="pearson")
cor_mat <- cor_X$r 
cor_p_mat <- cor_X$P
cor_p_mat[upper.tri(cor_p_mat, diag = FALSE)] <- p.adjust(cor_p_mat[upper.tri(cor_p_mat, diag = FALSE)], method = "BH")
cor_p_mat[lower.tri(cor_p_mat, diag = FALSE)] <- p.adjust(cor_p_mat[lower.tri(cor_p_mat, diag = FALSE)], method = "BH")
diag(cor_p_mat) <- 1
cor_p_mat[1:5,1:5]

extract_variable <- function(var_name){
  var_p <- adj_p_val[var_name,,drop=FALSE]
  IDs <- which(var_p<0.05)
  var_p <- var_p[1,IDs,drop=FALSE]
  var_r <- par_cor[var_name,,drop=FALSE]
  var_r <- var_r[1,IDs,drop=FALSE]
  res <- cbind(t(var_p), t(var_r))
  colnames(res) <- paste(colnames(res)[1], c("adj. p","par.cor."))
  return(res)
}
extract_variable("FOXA3_expr")
var_name <- "FOXA3_expr"
  
adj_p_val_expr_methylation <- adj_p_val[grep("_expr",rownames(adj_p_val)),grep("_methylation",rownames(adj_p_val))]
num_edges_expr_methylation <- apply(adj_p_val_expr_methylation, 1, function(x){return(sum(x<0.05))})
num_edges_expr_methylation <- sort(num_edges_expr_methylation, decreasing = TRUE)
head(num_edges_expr_methylation)

adj_p_val_expr_methylation_pearson <- cor_p_mat[grep("_expr",rownames(cor_p_mat)),grep("_methylation",rownames(cor_p_mat))]
cor_expr_methylation_pearson <- cor_mat[grep("_expr",rownames(cor_p_mat)),grep("_methylation",rownames(cor_p_mat))]
num_edges_expr_methylation_pearson <- apply(adj_p_val_expr_methylation_pearson, 1, function(x){return(sum(x<0.05))})
num_edges_expr_methylation_pearson <- sort(num_edges_expr_methylation_pearson, decreasing = TRUE)
head(num_edges_expr_methylation_pearson)
adj_p_val_expr_methylation_pearson[1:5,1:5]
dim(adj_p_val_expr_methylation_pearson)
#### analyze cis and trans associations ####
adj_p_val_expr_methylation_temp <- as.matrix(adj_p_val_expr_methylation)
sum(adj_p_val_expr_methylation_temp<0.05)
rownames(adj_p_val_expr_methylation_temp) <- gsub("_expr","",rownames(adj_p_val_expr_methylation_temp))
colnames(adj_p_val_expr_methylation_temp) <- gsub("_methylation","",colnames(adj_p_val_expr_methylation_temp))
features <- sort(intersect(rownames(adj_p_val_expr_methylation_temp), colnames(adj_p_val_expr_methylation_temp)))
length(features)
adj_p_val_expr_methylation_temp <- adj_p_val_expr_methylation_temp[features,features]  
sum(adj_p_val_expr_methylation_temp<0.05)
sum(diag(adj_p_val_expr_methylation_temp<0.05))
333/769
par_cor_temp <- as.matrix(par_cor)
rownames(par_cor_temp) <- gsub("_expr","",rownames(par_cor_temp))
colnames(par_cor_temp) <- gsub("_methylation","",colnames(par_cor_temp))
par_cor_temp <- par_cor_temp[features,features]
sum(diag(par_cor_temp)[diag(adj_p_val_expr_methylation_temp<0.05)]>0)
sum(diag(par_cor_temp)[diag(adj_p_val_expr_methylation_temp<0.05)]<0)
sum(diag(par_cor_temp)[diag(adj_p_val_expr_methylation_temp<0.05)]>0)/333
sum(diag(par_cor_temp)[diag(adj_p_val_expr_methylation_temp<0.05)]<0)/333


#### analyze cis and trans associations ####
adj_p_val_expr_methylation_pearson_temp <- as.matrix(adj_p_val_expr_methylation_pearson)
adj_p_val_expr_methylation_pearson_temp[1:5,1:5]
sum(adj_p_val_expr_methylation_pearson<0.05)
rownames(adj_p_val_expr_methylation_pearson_temp) <- gsub("_expr","",rownames(adj_p_val_expr_methylation_pearson_temp))
colnames(adj_p_val_expr_methylation_pearson_temp) <- gsub("_methylation","",colnames(adj_p_val_expr_methylation_pearson_temp))
features <- sort(intersect(rownames(adj_p_val_expr_methylation_pearson_temp), colnames(adj_p_val_expr_methylation_pearson_temp)))
length(features)
adj_p_val_expr_methylation_pearson_temp <- adj_p_val_expr_methylation_pearson_temp[features,features]  
sum(adj_p_val_expr_methylation_pearson_temp<0.05)
sum(diag(adj_p_val_expr_methylation_pearson_temp<0.05))
sum(diag(adj_p_val_expr_methylation_pearson_temp<0.05))/sum(adj_p_val_expr_methylation_pearson_temp<0.05)

cor_pearson_temp <- as.matrix(cor_mat)
rownames(cor_pearson_temp) <- gsub("_expr","",rownames(cor_pearson_temp))
colnames(cor_pearson_temp) <- gsub("_methylation","",colnames(cor_pearson_temp))
cor_pearson_temp <- cor_pearson_temp[features,features]
sum(diag(cor_pearson_temp)[diag(adj_p_val_expr_methylation_pearson_temp<0.05)]>0)/sum(diag(adj_p_val_expr_methylation_pearson_temp<0.05))
sum(diag(cor_pearson_temp)[diag(adj_p_val_expr_methylation_pearson_temp<0.05)]<0)/sum(diag(adj_p_val_expr_methylation_pearson_temp<0.05))

#### analyze cis and trans associations ####

threshold <- sort(abs(as.vector(cor_expr_methylation_pearson)), decreasing=TRUE)[769]
cor_expr_methylation_pearson_temp <- as.matrix(cor_expr_methylation_pearson)
rownames(cor_expr_methylation_pearson_temp) <- gsub("_expr","",rownames(cor_expr_methylation_pearson_temp))
colnames(cor_expr_methylation_pearson_temp) <- gsub("_methylation","",colnames(cor_expr_methylation_pearson_temp))
features <- sort(intersect(rownames(cor_expr_methylation_pearson_temp), colnames(cor_expr_methylation_pearson_temp)))
length(features)
cor_expr_methylation_pearson_temp <- cor_expr_methylation_pearson_temp[features,features]  
sum(abs(cor_expr_methylation_pearson_temp)>threshold)
sum(diag(abs(cor_expr_methylation_pearson_temp)>threshold))
sum(diag(abs(cor_expr_methylation_pearson_temp)>threshold))/sum(cor_expr_methylation_pearson_temp>threshold)

cor_pearson_temp <- as.matrix(cor_mat)
rownames(cor_pearson_temp) <- gsub("_expr","",rownames(cor_pearson_temp))
colnames(cor_pearson_temp) <- gsub("_methylation","",colnames(cor_pearson_temp))
cor_pearson_temp <- cor_pearson_temp[features,features]
sum(diag(cor_pearson_temp)[diag(adj_p_val_expr_methylation_pearson_temp<0.05)]>0)/sum(diag(adj_p_val_expr_methylation_pearson_temp<0.05))
sum(diag(cor_pearson_temp)[diag(adj_p_val_expr_methylation_pearson_temp<0.05)]<0)/sum(diag(adj_p_val_expr_methylation_pearson_temp<0.05))





res <- extract_variable("ZFP57_expr")#ZFP57 acts by controlling DNA methylation during the earliest multicellular stages of development at multiple imprinting control regions (ICRs) (PubMed:18622393, PubMed:30602440).
res <- res[order(abs(res[,2])),]

res <- extract_variable("ZNF334_expr")#ZNF334 hypermethylation in HCC survival relevant https://doi.org/10.1038/s41419-022-04895-6
#https://doi.org/10.1007/s00018-022-04295-1, 
res <- res[order(abs(res[,2])),]
res
gene_names_query <- gsub("_methylation", "", gsub("_expr","", rownames(res)))

res <- extract_variable("NR6A1_expr") #also known as GCNF, which interacts with DNMT3B,	DNA (cytosine-5)-methyltransferase 3B; Required for genome-wide de novo methylation and is essential for the establishment (STRING), https://doi.org/10.1016/j.bbrc.2006.04.007
res <- res[order(abs(res[,2])),]
res
gene_names_query <- gsub("_methylation", "", gsub("_expr","", rownames(res)))
adj_p_val[c(651, 2110),c(651, 2110)]
par_cor[c(651, 2110),c(651, 2110)]

res <- extract_variable("MYRFL_expr") 
res <- res[order(abs(res[,2])),]
res
gene_names_query <- gsub("_methylation", "", gsub("_expr","", rownames(res)))

getBM(attributes=c('hgnc_symbol','chromosome_name', 'start_position', 'end_position', 'strand'),
      filters=c('hgnc_symbol'),
      values=list(gene_names_query),
      mart=ensembl)


summarize_target_list <- function(target_list){
  summary_list = NULL
  for(i in target_list){
    temp_targets <- extract_variable(i)
    temp_targets <- temp_targets[-grep("_expr",rownames(temp_targets)),]
    rownames(temp_targets) <- gsub("_methylation", "", rownames(temp_targets))
    temp_targets <- temp_targets[order(abs(temp_targets[,2]), decreasing=TRUE),]
    temp_targets <-  cbind(temp_targets, sign_rho=paste(rownames(temp_targets), 
                                                        c("(-)","(+)")[as.integer(temp_targets[,2]>0)+1], sep=""))
    #print(temp_targets)
    
    summary_list <- rbind(summary_list, c(gsub("_expr","",i), nrow(temp_targets), paste(temp_targets[,"sign_rho"], collapse = ', ')))
  }
  colnames(summary_list) <- c("TF RNA", "#related methylation sites", "related methylation sites (+/-=edge sign)")
  return(summary_list)
} 

library(xtable)
xtable(summarize_target_list(c("ZFP57_expr", "ZNF334_expr", "NR6A1_expr", "MYRFL_expr")))

paste(c("A","B"), c("-","+"))
