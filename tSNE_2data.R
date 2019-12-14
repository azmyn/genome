library(data.table)
library(tictoc)
library(Rtsne)
library(NLP)
library(tm)
library(ggplot2)
library(ggbiplot)

# データ読み込み -----------------------------------------------------------------


# args <- commandArgs(trailingOnly = TRUE) #引数受け取り
tic()
data1 = fread(
  "/Users/azumi/Genome/luad.mut_1Mb_count.txt",
  stringsAsFactors = FALSE,
  encoding = "UTF-8"
)

data2 = fread(
  "/Users/azumi/Genome/lihc.mut_1Mb_count.txt",
  stringsAsFactors = FALSE,
  encoding = "UTF-8"
)


# mc3 = fread(args[1], stringsAsFactors = FALSE, encoding = "UTF-8")
# clinical = fread(
#   args[2],
#   stringsAsFactors = FALSE,
#   encoding = "UTF-8",
#   skip = 4
# )
toc()

print("データ読み込み完了")


# 行列作成 1 --------------------------------------------------------------------

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


mat = data.frame(matrix(rep(0, sum(chrmlen_1Mb)), nrow = 1))[numeric(0),]
ID =  data1$PATIENT_ID[1]
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

for (i in 1:nrow(data1)) {
  newID = data1$PATIENT_ID[i]
  if (newID != ID) {
    j = j + 1
    ID = newID
  }
  c = calc(data1$CHROMOSOME[i], data1$POSITION[i])
  mat[j, c] = data1$AMOUNT_OF_SUBSTITUTIONS[i]
}
# rownames(mat) = c(1:length(mat))

mat[is.na(mat)] <- 0
toc()

mat1 = mat

# 行列作成2 -------------------------------------------------------------------


mat = data.frame(matrix(rep(0, sum(chrmlen_1Mb)), nrow = 1))[numeric(0),]
ID =  data2$PATIENT_ID[1]
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

for (i in 1:nrow(data2)) {
  newID = data2$PATIENT_ID[i]
  if (newID != ID) {
    j = j + 1
    ID = newID
  }
  c = calc(data2$CHROMOSOME[i], data2$POSITION[i])
  mat[j, c] = data2$AMOUNT_OF_SUBSTITUTIONS[i]
}
# rownames(mat) = c(1:length(mat))

mat[is.na(mat)] <- 0
toc()

mat2 = mat

# 単語行列 ------------------------------------------------------------------

mat = rbindlist(list(mat1, mat2))
# mat2 = cbind(mat,as.data.table(c(rep("LUAD", nrow(mat1)),rep("LUSC",nrow(mat2)))) )

# tSNE --------------------------------------------------------------------
D = as.matrix(mat)
tic()
tsne = Rtsne(
  D,
  check_duplicates = FALSE,
  verbose = TRUE,
  initial_dims = 3126
)
toc()

tic()
png("~/Genome/luad_and_lihc.1Mb_tsne.png",
    width = 2000,
    height = 2000,
    )
colors = c(rep("2", nrow(mat1)),rep("3",nrow(mat2)))

plot(tsne$Y, t = 'n', main = "Rtsne")
text(tsne$Y, labels = as.character(as.factor(row.names(mat))),col = colors)
dev.off()
toc()


# PCA ---------------------------------------------------------------------
df = as.data.frame(mat)
df_f <- df[,apply(df, 2, var, na.rm=TRUE) != 0]

# biplot(x=dpca)
dpca = prcomp(df_f)
png("~/Genome/luad_and_lihc.1Mb_PCA.png",
    width = 2000,
    height = 2000)
print(
  ggbiplot(
    dpca,
    obs.scale = 1,
    var.scale = 1,
    groups = colors,
    ellipse = TRUE,
    circle = TRUE
  ) + ggtitle("LUAD(red) and LIHC(green)")
)
dev.off()
