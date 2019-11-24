library(data.table)


# データ読み込み -----------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE) #引数受け取り 1つ目がmutation, 2つ目がclinical
mut_dat = fread(args[1],stringsAsFactors=FALSE, encoding="UTF-8")
clin_dat = fread(args[2],stringsAsFactors=FALSE, encoding="UTF-8")

