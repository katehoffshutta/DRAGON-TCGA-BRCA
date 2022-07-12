library(aws.s3)
library(GenomicDataCommons)

gdc_set_cache(directory = "../../data/external/tcga_BRCA_gene_expression")

ge_manifest = files() %>%
    GenomicDataCommons::filter( cases.project.project_id == 'TCGA-BRCA') %>% 
    GenomicDataCommons::filter( type == 'gene_expression' ) %>%
    GenomicDataCommons::filter( data_type == 'Gene Expression Quantification') %>%
    GenomicDataCommons::filter( cases.demographic.gender == 'female') %>%
    manifest() %>%
    dplyr::filter(!grepl("splice_junctions",filename))

head(ge_manifest)
write.table(ge_manifest,file = "../../data/external/TCGA_BRCA_gene_expression_manifest.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Pulled files on 20220709
for(i in 1:length(ge_manifest$id))
{     
      options(warn=2)
      print(paste("Processing file:",i))
      
      fullpath = gdcdata(ge_manifest$id[[i]])
      dirname = names(fullpath)
      filename = tail(strsplit(fullpath,"/")[[1]],n=1)
      print(paste("Filename:",filename))
      
      put_object(file=fullpath,
      object = filename,
      bucket = paste(c("netzoo/supData/dragon/dragonDataPreprocessing/gene_expression/",dirname),collapse=""), 
      region="us-east-2",
      multipart=F)
      	
	# remove file from system
	splitpath = strsplit(fullpath,"/")[[1]]
	folderpath = paste(splitpath[-length(splitpath)],collapse="/")
	system2(command="rm",args=c("-r",folderpath))
}




