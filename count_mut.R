library(data.table)
library(tictoc)

tic()
mc3 = fread("/Users/azumi/Genome/mc3/mc3.v0.2.8.PUBLIC.maf.gz",stringsAsFactors=FALSE, encoding="UTF-8")
toc()