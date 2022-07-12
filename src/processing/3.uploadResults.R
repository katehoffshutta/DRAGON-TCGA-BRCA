library(aws.s3)
Sys.setenv("AWS_PROFILE" = "MFA")

put_object(file="../../data/processed/Pooled/dragon_mat.tsv",
      bucket = "netzoo/supData/dragon/dragonOutputFiles/Pooled",
      region="us-east-2",
      multipart=F)

put_object(file="../../data/processed/Pooled/dragon_adj_p.tsv",
      bucket = "netzoo/supData/dragon/dragonOutputFiles/Pooled",
      region="us-east-2",
      multipart=F)

put_object(file="../../data/processed/Pooled/dragon_raw_p.tsv",
      bucket = "netzoo/supData/dragon/dragonOutputFiles/Pooled",
      region="us-east-2",
      multipart=F)

