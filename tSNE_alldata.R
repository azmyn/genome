library(data.table)
library(tictoc)
library(Rtsne)
library(NLP)
library(tm)
library(ggplot2)
library(ggbiplot)
library(RColorBrewer) #カラーパレット
library(dplyr) #カラーパレット


mat = fread(
  "matrix.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)
mat = fread(
  "matrix_6types.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)
labels = fread(
  "matrix_labels.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)

label = labels[,1]
type = labels[,2]

# 行列作成(確率バージョン) -----------------------------------------------------------

sum = apply(mat, 1, sum)
prob <- function(x) {
  return(x / sum(x))
}
mat2 = as.matrix (mat)
mat_prob = matrix(rep(0, sum(chrmlen_1Mb) * nrow(mat2)), nrow = nrow(mat2))
for (i in 1:nrow(mat)) {
  for (j in 1:ncol(mat)) {
    mat_prob[i, j] = mat2[i, j] / sum[i]
  }
}
colnames(mat_prob) = cname

# tSNE --------------------------------------------------------------------
D = as.matrix(mat)
tic()
tsne = Rtsne(
  D,
  check_duplicates = FALSE,
  verbose = TRUE,
  initial_dims = nrow(mat)
)
toc()


tic()
png("~/Genome/tsne_alldata_6type.png",
    width = 2000,
    height = 2000,
    )

plot(tsne$Y, t = 'n', main = "Rtsne")
legend("bottomleft", legend = sort(unique(type)), col = c(1:25),pch = c(1:25)%/%6+15)
# text(tsne$Y, labels = as.character(as.factor(label)),col = type)
points(tsne$Y,col = type_num, pch = (type_num%/%6)+15)

dev.off()
toc()

# tSNE prob--------------------------------------------------------------------
D = as.matrix(mat_prob)
tic()
tsne = Rtsne(
  D,
  check_duplicates = FALSE,
  verbose = TRUE,
  initial_dims = nrow(mat_prob)
)
toc()
#プロット関連
toPoint <- function(factors) { 
  mapping <- c ("Biliary-AdenoCA"= 1,
                "Bone-Cart"= 2,
                "Bone-Epith"= 3,
                "Bone-Osteosarc"= 4,
                "Breast-AdenoCa"= 5,
                "Breast-DCIS"= 6,
                "Breast-LobularCa"= 7,
                "CNS-Medullo"= 8,
                "CNS-PiloAstro"= 9,
                "Eso-AdenoCa"= 10,
                "Head-SCC"= 11,
                "Kidney-RCC"= 12,
                "Liver-HCC"= 13,
                "Lymph-BNHL"= 14,
                "Lymph-CLL"= 15,
                "Lymph-NOS"= 16,
                "Myeloid-AML"= 17,
                "Myeloid-MDS"= 18,
                "Myeloid-MPN"= 19,
                "Ovary-AdenoCA"= 20,
                "Panc-AdenoCA"= 21,
                "Panc-Endocrine"= 22,
                "Prost-AdenoCA"= 23,
                "Skin-Melanoma"= 24,
                "Stomach-AdenoCA"= 25)
  mapping[as.character(factors)]
}
type_num <- unlist(lapply(type,toPoint))
k = unique(type)
s = c()
for (i in 1:25) {
  s = unlist(c(s,k[i]))
}
s=sort(s)
colPal3 <- colorRampPalette(brewer.pal(11, "Spectral"))
cp3 = colPal3(25)
#おしまい

png("~/Genome/tsne_alldata_6type_prob.png",
    width = 2000,
    height = 2000,
)
plot(tsne$Y, t = 'n', main = "Rtsne 頻度行列")
# legend("bottomleft", legend = s , col = cp3 ,pch = c(1:25)%/%6+15)
# # text(tsne$Y, labels = as.character(as.factor(label)),col = type)
# points(tsne$Y,col = cp3[type_num], pch = (type_num%/%6)+15 )

legend("bottomleft", legend = sort(unique(type)), col = c(1:25),pch = c(1:25)%/%6+15)
# text(tsne$Y, labels = as.character(as.factor(label)),col = type)
points(tsne$Y,col = type_num, pch = (type_num%/%6)+15)

dev.off()
