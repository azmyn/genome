library(data.table)
library(tictoc)
library(Rtsne)
library(NLP)
library(tm)



# データ読み込み -----------------------------------------------------------------


# args <- commandArgs(trailingOnly = TRUE) #引数受け取り
tic()
data = fread(
  "/Users/azumi/Dropbox/KU/shimolab_2019/genome/luad.mut_1Mb_count.txt",
  stringsAsFactors = FALSE,
  encoding = "UTF-8"
)
# clinical = fread(
#   "/Users/azumi/Genome/luad_tcga/data_bcr_clinical_data_patient.txt",
#   stringsAsFactors = FALSE,
#   encoding = "UTF-8",
#   skip = 4
# )

# mc3 = fread(args[1], stringsAsFactors = FALSE, encoding = "UTF-8")
# clinical = fread(
#   args[2],
#   stringsAsFactors = FALSE,
#   encoding = "UTF-8",
#   skip = 4
# )
toc()

print("データ読み込み完了")


# 行列作成 --------------------------------------------------------------------

chrmlen = c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415)
chrmlen_1Mb = chrmlen%/%1000000+1
sum(chrmlen_1Mb)

table(data$POSITION)


# 単語行列 ------------------------------------------------------------------

tm_corpus <- Corpus(VectorSource( ))
tdm <- TermDocumentMatrix(tm_corpus)
D = as.matrix(tdm)
setname <- function(x, y, z) {
  return (sprintf("%s_%s_%s", x, y, z))
}
col.name = mapply(setname, data$PATIENT_ID, data$CHROMOSOME, data$POSITION)

colnames(D) = as.vector(col.name)


# tSNE --------------------------------------------------------------------
tic()
tsne2 = Rtsne(t(D),
              check_duplicates = FALSE,
              verbose = TRUE,
              initial_dims = nrow(D)
)
toc()

tic()
tsne = Rtsne(t(D),
             check_duplicates = FALSE,
             verbose = TRUE)
toc()

tic()
png("~/Genome/1Mb_tsne.png",
    width = 1200,
    height = 1200)
plot(tsne$Y, t = 'n', main = "Rtsne")
text(tsne$Y, labels=as.character(as.factor(col.name)))
toc()

sum(D$AMOUNT_OF_SUBSTITUTIONS)