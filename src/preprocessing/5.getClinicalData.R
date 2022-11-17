library(dplyr)
library(GenomicDataCommons)
case_ids = cases() %>% GenomicDataCommons::filter(~ project.project_id == "TCGA-BRCA") %>% 
  GenomicDataCommons::filter(demographic.gender == 'female') %>%
  results_all() %>% ids()
clindat = gdc_clinical(case_ids)

demographics_small = clindat$demographic %>% 
  dplyr::select(c("race","gender","ethnicity","case_id"))
labeled_demographic = left_join(demographics_small,clindat$main,by="case_id")
write.table(labeled_demographic,
            file = "data/external/TCGA_BRCA_clinical.tsv",
            sep = "\t",
            row.names = F,
            quote = F)
