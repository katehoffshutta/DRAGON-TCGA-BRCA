library(GenomicDataCommons)

gdc_set_cache(directory = "data/external/tcga_BRCA_methylation")

ge_manifest = files() %>%
  GenomicDataCommons::filter( cases.project.project_id == 'TCGA-BRCA') %>%
  GenomicDataCommons::filter( type == 'methylation_beta_value' ) %>%
  GenomicDataCommons::filter( platform == "illumina human methylation 450") %>%
  GenomicDataCommons::filter( cases.demographic.gender == 'female') %>%
  manifest()
dim(ge_manifest)
# [1] 885   5

head(ge_manifest)
write.table(ge_manifest,file = "data/external/TCGA_BRCA_methylation_manifest.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Data in manuscript pulled 20220709
for(i in 300:length(ge_manifest$id))
{     
  options(warn=2)
  print(paste("Processing file:",i))
  
  fullpath = gdcdata(ge_manifest$id[[i]])
  dirname = names(fullpath)
  filename = tail(strsplit(fullpath,"/")[[1]],n=1)
  print(paste("Filename:",filename))
  
}


