library(ggplot2)
library(ggpubr)
library(gridExtra)
library(data.table)
library(dplyr)

# load the subtype data
meth = data.table(fread("data/interim/analysis_dataset.tsv",sep="\t",header=T))

comm5_genes = read.table("data/interim/comm5.tsv",sep="\t",header=T)


# more numbers for paper
length(grep("methylation",comm5_genes$gene))
# 25
length(grep("expr",comm5_genes$gene))
# 14

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

if(!dir.exists("reports/figures/communityDetection/community5"))
  dir.create("reports/figures/communityDetection/community5",recursive = T)

jpeg("reports/figures/communityDetection/community5/TFAP2A_methylation.jpeg",width=6,height=6,units="in",res=300)
methBoxplots("TFAP2A_methylation")
dev.off()

jpeg("reports/figures/communityDetection/community5/TFAP2B_methylation.jpeg",width=6,height=6,units="in",res=300)
methBoxplots("TFAP2B_methylation")
dev.off()

jpeg("reports/figures/communityDetection/community5/TFAP2C_methylation.jpeg",width=6,height=6,units="in",res=300)
methBoxplots("TFAP2C_methylation")
dev.off()

jpeg("reports/figures/communityDetection/community5/TFAP2A_expression.jpeg",width=6,height=6,units="in",res=300)
methBoxplots("TFAP2A_expr")
dev.off()

jpeg("reports/figures/communityDetection/community5/TFAP2B_expression.jpeg",width=6,height=6,units="in",res=300)
methBoxplots("TFAP2B_expr")
dev.off()

jpeg("reports/figures/communityDetection/community5/TFAP2C_expression.jpeg",width=6,height=6,units="in",res=300)
methBoxplots("TFAP2C_expr")
dev.off()

kruskal.test(TFAP2A_methylation ~ Subtype_mRNA ,data = meth)
# Kruskal-Wallis chi-squared = 20.863, df = 4, p-value = 0.0003371

kruskal.test(TFAP2B_methylation ~ Subtype_mRNA ,data = meth)
# Kruskal-Wallis chi-squared = 152.45, df = 4, p-value < 2.2e-16

kruskal.test(TFAP2C_methylation ~ Subtype_mRNA ,data = meth)
# Kruskal-Wallis chi-squared = 6.3547, df = 4, p-value = 0.1742

kruskal.test(TFAP2A_expr ~ Subtype_mRNA ,data = meth)
# Kruskal-Wallis chi-squared = 60.902, df = 4, p-value = 1.875e-12

kruskal.test(TFAP2B_expr~ Subtype_mRNA ,data = meth)
# Kruskal-Wallis chi-squared = 175.43, df = 4, p-value < 2.2e-16

kruskal.test(TFAP2C_expr~ Subtype_mRNA ,data = meth)
# Kruskal-Wallis chi-squared = 51.362, df = 4, p-value = 1.876e-10

