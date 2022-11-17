library(GenomicDataCommons)

gdc_set_cache(directory = "data/external/tcga_BRCA_gene_expression")

ge_manifest = files() %>%
    GenomicDataCommons::filter( cases.project.project_id == 'TCGA-BRCA') %>% 
    GenomicDataCommons::filter( type == 'gene_expression' ) %>%
    GenomicDataCommons::filter( data_type == 'Gene Expression Quantification') %>%
    GenomicDataCommons::filter( cases.demographic.gender == 'female') %>%
    manifest() %>%
    dplyr::filter(!grepl("splice_junctions",filename))

write.table(ge_manifest,file = "data/external/TCGA_BRCA_gene_expression_manifest.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Pulled manuscript files on 20220709
for(i in 1:length(ge_manifest$id))
{     
      options(warn=2)
      print(paste("Processing file:",i))
      
      fullpath = gdcdata(ge_manifest$id[[i]])
      dirname = names(fullpath)
      filename = tail(strsplit(fullpath,"/")[[1]],n=1)
      print(paste("Filename:",filename))

}




