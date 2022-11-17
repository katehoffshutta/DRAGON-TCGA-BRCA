# 20220420
# This script performs thresholding based on FDR, community detection, and gsea
# The input is (1) a dragon adjacency matrix (2) a dragon adj p matrix 
# (3) a pathway list (gmt format)
# The output is (0) diagnostic plots, including FDR dist. and dist of edge weights
# (1) a hairball graph (2) a table consisting of all the gsea results
# with FDRp < 0.05 for every community (3) N graphs of individual communities, where N is the number
# of communities that have at least one pathway 

# libraries
library(data.table)
library(dplyr)
library(fgsea)
library(igraph)
library(tidyverse)

prefix = "Pooled5"
thres = 0.005 # arbitrary

res = data.frame(fread(paste(c("data/processed",prefix,"dragon_mat.tsv"),collapse="/"),
                       sep="\t",header=T), row.names = 1)
adj_p = data.frame(fread(paste(c("data/processed",prefix,"dragon_adj_p.tsv"),collapse="/"),
                         sep="\t",header=T),row.names = 1)
input = data.frame(fread(paste(c("data/processed",prefix,"dragon_input_mat.tsv"),collapse="/"),
                         sep="\t",header=T),row.names = 1)

# get numbers for paper
nExp = length(grep("_expr",names(input)))
nMeth = length(grep("_methylation",names(input)))
dim(input)
print(nExp)
#1311
print(nMeth)
#1557
print(nMeth + nExp)
#2868

# for how many genes do we have both?
exprGenesOnly = sapply(names(input)[grep("expr",names(input))],function(x){strsplit(x,"_")[[1]][1]})
methGenesOnly = sapply(names(input)[grep("methylation",names(input))],function(x){strsplit(x,"_")[[1]][1]})

length(intersect(exprGenesOnly,methGenesOnly))
#1280 in both

meth_not_exp = setdiff(methGenesOnly,exprGenesOnly)
exp_not_meth = setdiff(exprGenesOnly,methGenesOnly)
all_diff = c(meth_not_exp, exp_not_meth)
write.table(all_diff,
            "data/interim/gene_names_diff.tsv",
            sep="\t",row.names=F,col.names = "geneName")

folder = paste("reports/figures",prefix,sep="/")
if (!file.exists(folder))
  dir.create(folder,recursive = T)

pdf(paste(c("reports/figures",prefix,"diagnostics.pdf"),collapse="/"))
hist(as.matrix(adj_p),main="DRAGON edge FDR distribution", xlab="FDR")
hist(as.matrix(adj_p[adj_p < 0.05]),main="DRAGON edge FDR distribution \n FDR < 0.05", xlab="FDR")
hist(as.matrix(res),breaks=100, xlab="partial correlation", 
               main="distribution of partial correlations")
dev.off()
sum(adj_p < 0.005)

jpeg("reports/figures/FDRdist.jpeg",width=8,height=6,units="in",res=300)
par(mfrow=c(1,2))
hist(as.matrix(adj_p),main="DRAGON edge FDR distribution", xlab="FDR")
hist(as.matrix(adj_p[adj_p < 0.05]),main="DRAGON edge FDR distribution \n FDR < 0.05", xlab="FDR")
dev.off()
# threshold based on input FDRp
res_thres = res
res_thres[adj_p > thres] = 0

# sanity check: these two numbers should be the same
sum(res_thres != 0)
sum(adj_p < 0.005) - nMeth - nExp

# make hairball graph
mygraph = graph_from_adjacency_matrix(as.matrix(res_thres),mode="undirected",diag=F,weighted=T)
mysmallgraph = induced_subgraph(mygraph,degree(mygraph)>0)
label_colors = rep(NA,length(V(mysmallgraph)))
label_colors[grep("meth",V(mysmallgraph)$name)] = "darkorange"
label_colors[grep("expr",V(mysmallgraph)$name)] = "turquoise4"
edge_colors = ifelse(E(mysmallgraph)$weight > 0,"red","blue")

set.seed(2022)
mylayout = layout_with_graphopt(mysmallgraph)

jpeg(paste(c("reports/figures",prefix,"hairball.jpeg"),collapse="/"),width=10,height=10,units="in",res=300)

plot(mysmallgraph,
     vertex.label=NA,
     edge.color = edge_colors,
     edge.width = 10*abs(E(mysmallgraph)$weight),
     vertex.size= 1.5,
     vertex.color = label_colors,
     vertex.frame.color= label_colors,
     vertex.label.cex=0.2,
     vertex.label.color=label_colors,
     vertex.label.dist=0.2,
     layout=mylayout,
     main= paste("DRAGON BRCA GGM: FDRp < 0.005",prefix,sep="\n"))
legend(-0.55,-1.1,c("methylation","expression"),col=c("darkorange","turquoise4"),pch=c(20,20), pt.cex = 3, cex = 1.2)
legend(0,-1.1,c("positive partial cor","negative partial cor"),col=c("red","blue"),lty=c(1,1),lwd=c(6,6), cex = 1.2)
dev.off()

# some stats for paper
length(E(mysmallgraph))
# 3631 edges
length(V(mygraph))
# 2868 total nodes
length(V(mysmallgraph))
# 2106 nodes with degree > 0
length(grep("methylation",V(mygraph)$name))
# 1557
length(grep("methylation",V(mysmallgraph)$name))
# 1168; 1168/1557=0.75
length(grep("expr",V(mygraph)$name))
# 1311
length(grep("expr",V(mysmallgraph)$name))
# 938; 938/1311 = 0.72

# do community detection
# Use cluster walktrap and find max modularity
bestMod=0
bestSteps=0

for(i in 1:50)
{
  set.seed(42)
  comms = cluster_walktrap(mysmallgraph, steps=i,
                         weights=abs(E(mysmallgraph)$weight))
  thisMod = modularity(mysmallgraph, comms$membership)
  if(thisMod > bestMod)
  {
    bestMod = thisMod
    bestSteps = i
  }
}

set.seed(42)
comms = cluster_walktrap(mysmallgraph,steps=bestSteps,
                         weights=abs(E(mysmallgraph)$weight))

# how does this compare to louvain
comms_louvain = cluster_louvain(mysmallgraph,
                         weights=abs(E(mysmallgraph)$weight))

comms_fg= cluster_fast_greedy(mysmallgraph,
                                weights=abs(E(mysmallgraph)$weight))

comms_eigen = cluster_leading_eigen(mysmallgraph,
                               weights=abs(E(mysmallgraph)$weight))

modularity(mysmallgraph,membership=comms$membership)
modularity(mysmallgraph,membership=comms_louvain$membership)
modularity(mysmallgraph,membership=comms_fg$membership)
modularity(mysmallgraph,membership=comms_eigen$membership)

# comms_fg is the best from modularity perspective
comms = comms_fg

# get some descriptions for the paper
length(unique(comms$membership))
# 169 communities
summary(factor(summary(factor(comms$membership),maxsum=300)))
#2   3   4   5   6   7   9  10  11  12  13  14 
#69  27  14   9   5   6   2   1   2   2   3   1 
#15  16  17  20  21  23  29  30  32  33  39  40 
#1   1   1   1   1   2   1   2   1   1   1   2 
#41  44  45  50  53  55  60  62  72  83  98 101 
#1   1   1   1   1   1   1   1   1   1   1   1 
#415 
#1 

# hairball comms plots
jpeg(paste(c("reports/figures",prefix,"hairball_comms.jpeg"),collapse="/"),width=10,height=10,units="in",res=300)
set.seed(2022)

# using code from https://stackoverflow.com/questions/16390221/how-to-make-grouped-layout-in-igraph

edge.weights <- function(community, network, weight.within = 100, weight.between = 1) {
  bridges <- igraph::crossing(communities = community, graph = network)
  weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
  return(weights) 
}

mysmallcommgraph = mysmallgraph
E(mysmallcommgraph)$weights = edge.weights(comms,mysmallgraph)
mysmallcommgraph = induced_subgraph(mysmallgraph, get_indx)
E(mysmallcommgraph)$weight = edge.weights(comms,mysmallcommgraph)
get_indx = which(comms$membership %in% c(6,14))
set.seed(46)
plot(mysmallcommgraph,
     vertex.label=NA,
     mark.groups = list(which(comms$membership[get_indx]==6),
                        which(comms$membership[get_indx]==14)),
     #edge.color = edge_colors[get_indx],
     #edge.width = 10*abs(E(mysmallcommgraph)$weight),
     vertex.size= 3,
     vertex.color = label_colors[get_indx],
     vertex.frame.color= label_colors[get_indx],
     vertex.label.cex=0.2,
     vertex.label.color=label_colors[get_indx],
     vertex.label.dist=0.2,
     layout=layout_with_fr,
     main="DRAGON Communities")
#legend("bottomleft",c("methylation","expression"),col=c("darkorange","turquoise4"),pch=c(20,20))
#legend("bottomright",c("positive partial cor","negative partial cor"),col=c("red","blue"),pch=c(20,20))
dev.off()

# get numbers for writeup
sum(summary(factor(comms$membership),maxsum = 300)==1)
sum(summary(factor(comms$membership),maxsum = 300)==2)
sum(summary(factor(comms$membership),maxsum = 300)==3)
sum(summary(factor(comms$membership),maxsum = 300)==4)
sum(summary(factor(comms$membership),maxsum = 300)>=5)


# do enrichment analysis
download.file("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/c2.cp.reactome.v7.5.1.symbols.gmt",
              destfile = "data/external/c2.cp.reactome.v7.5.1.symbols.gmt")

pathways_reactome = gmtPathways("data/external/c2.cp.reactome.v7.5.1.symbols.gmt") #../../data/external/c5.go.bp.v7.5.1.symbols.gmt")
pathways_methylation= lapply(pathways_reactome,function(x){paste(x,"methylation",sep="_")})
pathways_expression = lapply(pathways_reactome,function(x){paste(x,"expr",sep="_")})
pathways_both = lapply(pathways_reactome,function(x){a=paste(x,"methylation",sep="_"); b= paste(x,"expr",sep="_"); return(c(a,b))})

enrichmentResults = list()
enrichmentResultsMeth = list()
enrichmentResultsExpr = list()

for(i in 1:length(unique(comms$membership)))
{
  if(i%%10 == 0) print(paste("Testing community",i))
  enrichmentResults[[i]]=fora(pathways_both,
                              genes = comms$names[comms$membership == i],
                              universe = comms$names,
                              minSize = 3)
  enrichmentResultsMeth[[i]]=fora(pathways_methylation,
                              genes = comms$names[comms$membership == i],
                              universe = comms$names,
                              minSize = 3)
  enrichmentResultsExpr[[i]]=fora(pathways_expression,
                              genes = comms$names[comms$membership == i],
                              universe = comms$names,
                              minSize = 3)
}

# filter for communities enriched w/FDRp < 0.05 
fdrThres = 0.05

enrichedComms = which(sapply(enrichmentResults,function(x){sum(x$padj < fdrThres)})>0)
enrichedCommsMeth = which(sapply(enrichmentResultsMeth,function(x){sum(x$padj < fdrThres)})>0)
enrichedCommsExpr = which(sapply(enrichmentResultsExpr,function(x){sum(x$padj < fdrThres)})>0)

# consider any community in any of these is enriched

sigResults = enrichmentResults[[enrichedComms[[1]]]] %>% 
  dplyr::filter(padj < fdrThres)
sigResults$community = enrichedComms[[1]]
sigResults$commSize = length(comms$names[comms$membership == enrichedComms[[1]]])
sigResults$method = "both"

for(i in 2:length(enrichedComms))
{
  c = enrichedComms[[i]]
  theseResults = enrichmentResults[[c]] %>% dplyr::filter(padj < fdrThres)
  theseResults$community = c
  theseResults$commSize = length(comms$names[comms$membership == c])
  theseResults$method = "both"
  sigResults= rbind.data.frame(sigResults,theseResults)
}

for(i in 1:length(enrichedCommsMeth))
{
  c = enrichedCommsMeth[[i]]
  theseResults = enrichmentResultsMeth[[c]] %>% dplyr::filter(padj < fdrThres)
  theseResults$community = c
  theseResults$commSize = length(comms$names[comms$membership == c])
  theseResults$method = "methylation"
  sigResults= rbind.data.frame(sigResults,theseResults)
}

for(i in 1:length(enrichedCommsExpr))
{
  c = enrichedCommsExpr[[i]]
  theseResults = enrichmentResultsExpr[[c]] %>% dplyr::filter(padj < fdrThres)
  theseResults$community = c
  theseResults$commSize = length(comms$names[comms$membership == c])
  theseResults$method = "expression"
  sigResults= rbind.data.frame(sigResults,theseResults)
}

sigResults$overlapGenes = sapply(sigResults$overlapGenes,function(x){paste(x,collapse=";")})
# turn overlap into a string using collapse before applyingcharacter
write.table(apply(sigResults,2,as.character),paste(c("reports/figures/brcaCommsTable_",prefix,".tsv"),collapse=""),
            row.names=F,quote=F,sep="\t")

sigComms = sigResults %>% group_by(community,commSize,method) %>% 
  summarize(nsig = sum(padj < fdrThres)) %>% 
  dplyr::filter(commSize > 4) %>% dplyr::select(community,commSize,method) %>%
  dplyr::filter(method == "both")

# number for manuscript

length(unique(sigComms$community))
#10 communities are enriched at FDR < 0.05

pdf(paste(c("reports/figures/brcaComms_all59_",prefix,".pdf"),collapse=""),width=10,height=10)

for(i in 1:nrow(sigComms))
{
  theseNodes = comms$names[comms$membership == sigComms$community[i]]
  mysmallgraph = induced_subgraph(mygraph,V(mygraph)$name %in% theseNodes)
  label_colors = rep(NA,length(V(mysmallgraph)))
  label_colors[grep("meth",V(mysmallgraph)$name)] = "darkorange"
  label_colors[grep("expr",V(mysmallgraph)$name)] = "turquoise4"
  label_short = str_split(V(mysmallgraph)$name,"_",simplify=T)[,1]
  edge_colors = ifelse(E(mysmallgraph)$weight > 0,"red","blue")
  
  plot(mysmallgraph,
       edge.color = edge_colors,
       edge.width = 3*abs(E(mysmallgraph)$weight),
       vertex.size= 3,
       vertex.color = label_colors,
       vertex.frame.color= label_colors,
       vertex.label = label_short,
       vertex.label.cex=0.8,
       vertex.label.color=label_colors,
       vertex.label.dist=0.8,
       layout = layout_with_graphopt,
       main=paste(c("Community",sigComms$community[i],"Method",sigComms$method[i])))
  legend("bottomleft",c("methylation","expression"),col=c("darkorange","turquoise4"),pch=c(20,20))
  legend("bottomright",c("positive partial cor","negative partial cor"),col=c("red","blue"),pch=c(20,20))
}

dev.off()

# make individual plots for the special communities
# community of size  >4, overlap of size >3, enriched in >3 reactome pathways
# plus community 36 because it's its own little button 
# criteria for interesting: commSize > 4, 
# overlap size >3, 
# number of enriched paths > 3

interesting = sigResults %>% dplyr::filter(commSize > 4, overlap > 3) %>%
  group_by(community,commSize,method) %>% 
  summarize(nsig = sum(padj < 0.05)) %>% 
  dplyr::filter(nsig > 3) %>% dplyr::select(community,commSize,method)

# list by inspection
special = c(5,38,41)

for(s in special)
{
  print(sigResults %>% dplyr::filter(community == s) %>% 
    summarize(nsig = sum(padj < 0.05)))
}

# community 5: 12
# community 38: 25
# community 41: 14

sigResults %>% filter(community == 38) %>% 
  arrange(padj) %>% select(pathway) %>% unique() %>% nrow()

set.seed(34)

for(i in special)
{
  # write membership table for community 14 for analysis
  thisComm = data.frame("gene"=comms$names[comms$membership == i])
  write.table(thisComm,file=paste(c("data/interim/comm",i,".tsv"),collapse=""),row.names=F)
  
  filename = paste(c("reports/figures/brcaComms_enriched_",prefix,"_",i,".jpeg"),collapse="")
  jpeg(filename, width=10,height=10,units="in",res=300)
  theseNodes = comms$names[comms$membership == i]
  mysmallgraph = induced_subgraph(mygraph,V(mygraph)$name %in% theseNodes)
  label_colors = rep(NA,length(V(mysmallgraph)))
  label_colors[grep("meth",V(mysmallgraph)$name)] = "darkorange"
  label_colors[grep("expr",V(mysmallgraph)$name)] = "turquoise4"
  label_short = str_split(V(mysmallgraph)$name,"_",simplify=T)[,1]
  edge_colors = ifelse(E(mysmallgraph)$weight > 0,"red","blue")
  
  plot(mysmallgraph,
       edge.color = edge_colors,
       edge.width = 10*abs(E(mysmallgraph)$weight),
       vertex.size= 4,
       vertex.color = label_colors,
       vertex.frame.color= label_colors,
       vertex.label = label_short,
       vertex.label.font = 2,
       vertex.label.cex=1.1,
       vertex.label.color=label_colors,
       vertex.label.dist=1.1,
       layout = layout_with_graphopt)
  title(paste("Community",i), cex.main=3)
  legend(-1.25,-0.75,cex = 1.2,c("methylation","expression"),col=c("darkorange","turquoise4"),pch=c(20,20),pt.cex=c(3,3))
  legend(-1.25,-1,cex = 1.2,c("positive partial cor","negative partial cor"),col=c("red","blue"),lty=c(1,1),lwd=c(4,4))
  dev.off()
}

# more numbers for paper

length(comms$names[comms$membership==5])
# 39 
length(grep("_methylation",comms$names[comms$membership==5]))
# 25
length(grep("_expr",comms$names[comms$membership==5]))
# 14
sigResults %>% filter(community == 5 & overlap >=5) %>% select(overlapGenes)
# community 2
length(comms$names[comms$membership==38])
length(grep("_methylation",comms$names[comms$membership==38]))

# community 4
length(comms$names[comms$membership==41])
length(grep("_methylation",comms$names[comms$membership==41]))

# community 8 

# community 15

# community 39

# community 49

