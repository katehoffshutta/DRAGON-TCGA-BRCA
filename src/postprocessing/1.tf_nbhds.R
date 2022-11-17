# plot neighborhoods from expression of the following TFS
# ZFP57
# ZNF334
# NR6A1
# MYRFL
library(data.table)
library(dplyr)
library(igraph)
library(tidyverse)

prefix = "Pooled5"

res = data.frame(fread(paste(c("data/processed",prefix,"dragon_mat.tsv"),collapse="/"),
                       sep="\t",header=T), row.names = 1)
adj_p = data.frame(fread(paste(c("data/processed",prefix,"dragon_adj_p.tsv"),collapse="/"),
                         sep="\t",header=T),row.names = 1)

thres = 0.05
res_thres = res
res_thres[adj_p > thres] = 0

myGraph = graph_from_adjacency_matrix(as.matrix(res_thres),weighted=T,diag=F,mode="undirected")
V(myGraph)$name = str_replace(V(myGraph)$name,"_"," ")
V(myGraph)$name = str_replace(V(myGraph)$name,"expr","gene expression")

set.seed(79)
jpeg("reports/figures/ZFP57_hub.jpeg",width=10,height=10,units="in",res=300)
zfp57 = make_ego_graph(myGraph,order=1,nodes="ZFP57 gene expression")[[1]]
label_colors = rep(NA,length(V(zfp57)))
label_colors[grep("meth",V(zfp57)$name)] = "darkorange"
label_colors[grep("expr",V(zfp57)$name)] = "turquoise4"
plot(zfp57,
     edge.width= 20*abs(E(zfp57)$weight),
     edge.color = ifelse(E(zfp57)$weight > 0, "red","blue"),
     vertex.label.dist = c(rep(2,13), 5),
     vertex.label = V(zfp57)$name,
     vertex.label.cex = 1.2,
     vertex.label.degree=-pi/2,
     #vertex.label.font = 2,
     vertex.label.color = "black",
     vertex.size = 5+10*log(degree(zfp57),base=1.8),
     vertex.color=label_colors,
     layout=layout_with_gem)
title("ZFP57 Gene Expression Hub", cex.main=2)
dev.off()

jpeg("reports/figures/ZNF334_hub.jpeg",width=10,height=10,units="in",res=300)
znf334 = make_ego_graph(myGraph,order=1,nodes="ZNF334 gene expression")[[1]]
label_colors = rep(NA,length(V(znf334)))
label_colors[grep("meth",V(znf334)$name)] = "darkorange"
label_colors[grep("expr",V(znf334)$name)] = "turquoise4"
plot(znf334,
     edge.width= 20*abs(E(znf334)$weight),
     edge.color = ifelse(E(znf334)$weight > 0, "red","blue"),
     vertex.label.dist = c(rep(2,12), 5),
     vertex.label.cex = 1.2,
     vertex.label.degree=-pi/2,
     #vertex.label.font = 2,
     vertex.label.color = "black",
     vertex.size = 5+10*log(degree(znf334),base=1.8),
     vertex.color=label_colors,
     layout=layout_with_gem)
title("ZNF334 Gene Expression Hub",cex.main=2)
dev.off()

jpeg("reports/figures/NR6A1_hub.jpeg",width=10,height=10,units="in",res=300)
nr6a1= make_ego_graph(myGraph,order=1,nodes="NR6A1 gene expression")[[1]]
label_colors = rep(NA,length(V(nr6a1)))
label_colors[grep("meth",V(nr6a1)$name)] = "darkorange"
label_colors[grep("expr",V(nr6a1)$name)] = "turquoise4"
plot(nr6a1,
     edge.width= 20*abs(E(nr6a1)$weight),
     edge.color = ifelse(E(nr6a1)$weight > 0, "red","blue"),
     vertex.label.dist = c(rep(2,9),1.5,2,2,2,5),
     vertex.label.cex = 1.2,
     vertex.label.degree=-pi/2,
     #vertex.label.font = 2,
     vertex.label.color = "black",
     vertex.size = 5+10*log(degree(nr6a1),base=1.8),
     vertex.color=label_colors,
     layout=layout_with_gem)
title("NR6A1 Gene Expression Hub",cex.main=2)
dev.off()

jpeg("reports/figures/MYRFL_hub.jpeg",width=10,height=10,units="in",res=300)
myrfl= make_ego_graph(myGraph,order=1,nodes="MYRFL gene expression")[[1]]
label_colors = rep(NA,length(V(myrfl)))
label_colors[grep("meth",V(myrfl)$name)] = "darkorange"
label_colors[grep("expr",V(myrfl)$name)] = "turquoise4"
plot(myrfl,
     edge.width= 20*abs(E(myrfl)$weight),
     edge.color = ifelse(E(myrfl)$weight > 0, "red","blue"),
     vertex.label.dist = c(rep(2,9),2,6,2,2),
     vertex.label.cex = 1.2,
     vertex.label.degree=-pi/2,
     #vertex.label.font = 2,
     vertex.label.color = "black",
     vertex.size = 5+10*log(degree(myrfl),base=1.8),
     vertex.color=label_colors,
     layout=layout_with_gem)
title("MYRFL Gene Expression Hub",cex.main=2)

dev.off()