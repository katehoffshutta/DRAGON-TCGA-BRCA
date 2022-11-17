# look at differential methylation of ZNFs by subtypes
# run full-postprocess.R before this to get the comms object
# load the methylation data

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(data.table)
library(dplyr)

# load the subtype data
meth = data.table(fread("data/interim/analysis_dataset.tsv",sep="\t",header=T))

comm38_genes = read.table("data/interim/comm38.tsv",sep="\t",header=T)


# more numbers for paper
length(grep("methylation",comm38_genes$gene))
# 6
length(grep("expr",comm38_genes$gene))
# 34

my_comparisons <- list( c("Basal","Her2"),
                        c("Basal","LumA"),
                        c("Basal","LumB"),
                        c("Basal","Normal"),
                        c("LumA","Her2"),
                        c("LumA","LumB"),
                        c("LumA","Normal"),
                        c("LumB","Her2"),
                        c("LumB","Normal"),
                        c("Normal","Her2"))

methBoxplots = function(myGene)
{
  ggboxplot(meth, x = "Subtype_mRNA", y = myGene, order=c("Basal","Her2","LumA","LumB","Normal"))+
    stat_compare_means(comparisons=my_comparisons) + 
    stat_compare_means(vjust=-20,color="red")
}

if(!dir.exists("reports/figures/communityDetection/community38/"))
  dir.create("reports/figures/communityDetection/community38/")

jpeg("reports/figures/communityDetection/community38/KLF6_methylation.jpeg",width=6,height=6,units="in",res=300)
methBoxplots("KLF6_methylation")
dev.off()

jpeg("reports/figures/communityDetection/community38/NR4A2_methylation.jpeg",width=6,height=6,units="in",res=300)
methBoxplots("NR4A2_methylation")
dev.off()

jpeg("reports/figures/communityDetection/community38/ATF1_methylation.jpeg",width=6,height=6,units="in",res=300)
methBoxplots("ATF1_methylation")
dev.off()

kruskal.test(KLF6_methylation ~ Subtype_mRNA ,data = meth)
# Kruskal-Wallis chi-squared = 93.037, df = 4, p-value < 2.2e-16

kruskal.test(NR4A2_methylation ~ Subtype_mRNA ,data = meth)
# Kruskal-Wallis chi-squared = 17.984, df = 4, p-value = 0.001243

kruskal.test(ATF1_methylation ~ Subtype_mRNA ,data = meth)
# Kruskal-Wallis chi-squared = 6.6151, df = 4, p-value = 0.1577



