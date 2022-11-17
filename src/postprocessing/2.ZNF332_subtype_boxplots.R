# ZNF334 subtype boxplots
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(data.table)
library(dplyr)

# load the subtype data
meth = data.table(fread("data/interim/analysis_dataset.tsv",sep="\t",header=T))


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

jpeg("reports/figures/ZNF334_methylation.jpeg",width=6,height=6,units="in",res=300)
methBoxplots("ZNF334_methylation")
dev.off()
