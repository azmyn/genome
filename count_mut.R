library(data.table)
library(tictoc)

args <- commandArgs(trailingOnly=TRUE) #引数受け取り

tic()
# mc3 = fread("/Users/azumi/Genome/mc3/mc3.v0.2.8.PUBLIC.maf.gz",stringsAsFactors=FALSE, encoding="UTF-8")
mc3 = fread(args[1],stringsAsFactors=FALSE, encoding="UTF-8")
toc()

test = mc3[1:1000,]

print(test)

