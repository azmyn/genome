library(data.table)
library(tictoc)
library(Rtsne)
library(NLP)
library(tm)
library(ggplot2)
library(ggbiplot)

# データ読み込み -----------------------------------------------------------------


args <- commandArgs(trailingOnly = TRUE) #引数受け取り

# data = fread(
#   "/Users/azumi/Genome/brca.mut_1Mb_count.txt",
#   stringsAsFactors = FALSE,
#   encoding = "UTF-8"
# )
# clinical = fread(
#   "/Users/azumi/Genome/luad_tcga/data_bcr_clinical_data_patient.txt",
#   stringsAsFactors = FALSE,
#   encoding = "UTF-8",
#   skip = 4
# )
data = fread(args[1], stringsAsFactors = FALSE, encoding = "UTF-8")
print("データ読み込み完了")

# 行列作成 --------------------------------------------------------------------

chrmlen = c(
  248956422,
  242193529,
  198295559,
  190214555,
  181538259,
  170805979,
  159345973,
  145138636,
  138394717,
  133797422,
  135086622,
  133275309,
  114364328,
  107043718,
  101991189,
  90338345,
  83257441,
  80373285,
  58617616,
  64444167,
  46709983,
  50818468,
  156040895,
  57227415
)
chrmlen_1Mb = ceiling(chrmlen / 1000000 + 1)


mat = data.frame(matrix(rep(0, sum(chrmlen_1Mb)), nrow = 1))[numeric(0), ]
ID =  data$PATIENT_ID[1]
j = 1
calc = function(x, y) {
  if (x == 'X') {
    x = 23
  }
  if (x == 'Y') {
    x = 24
  }
  if (x == 1) {
    return(y + 1)
  } else{
    x = as.numeric(x)
    return (sum(chrmlen_1Mb[1:(x - 1)]) + y + 1)
  }
}
tic()

for (i in 1:nrow(data)) {
  newID = data$PATIENT_ID[i]
  if (newID != ID) {
    j = j + 1
    ID = newID
  }
  c = calc(data$CHROMOSOME[i], data$POSITION[i])
  mat[j, c] = data$AMOUNT_OF_SUBSTITUTIONS[i]
}
# rownames(mat) = c(1:length(mat))
# setname <- function(x, y) {
#   return (sprintf("%s_%s", x, y))
# }
# col.name = mapply(setname,c(1:22,'X','Y'), chrmlen_1Mb)
mat[is.na(mat)] <- 0
toc()

# 単語行列 ------------------------------------------------------------------

# tm_corpus <- Corpus(VectorSource(data$MUTATIONS))
# tdm <- TermDocumentMatrix(tm_corpus)
# D = as.matrix(tdm)
# setname <- function(x, y, z) {
#   return (sprintf("%s_%s_%s", x, y, z))
# }
# col.name = mapply(setname, data$PATIENT_ID, data$CHROMOSOME, data$POSITION)
# colnames(D) = as.vector(col.name)


# tSNE --------------------------------------------------------------------
D = as.matrix(mat)
tic()
tsne = Rtsne(
  D,
  check_duplicates = FALSE,
  verbose = TRUE,
  initial_dims = ncol(D)
)
toc()

file.name = sprintf("~/Genome/%s_1Mb_tsne.png", substr(args[1], 1, 4))
png(file.name,
    width = 2000,
    height = 2000)
plot(tsne$Y, t = 'n', main = sprintf("tSNE(%s)", substr(args[1], 1, 4)))
text(tsne$Y, labels = as.character(as.factor(row.names(D))))
dev.off()


# PCA ---------------------------------------------------------------------

dpca = prcomp(t(D))
file.name = sprintf("~/Genome/%s_1Mb_PCA.png", substr(args[1], 1, 4))
png(file.name,
    width = 2000,
    height = 2000)
print(
  ggbiplot(
    dpca,
    obs.scale = 1,
    var.scale = 1,
    # groups = dat$SMOKING_HISTORY,
    ellipse = TRUE,
    circle = TRUE
  ) + ggtitle(sprintf("PCA(%s)", substr(args[1], 1, 4)))
)
dev.off()
