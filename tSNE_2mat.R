library(data.table)
library(tictoc)
library(umap)
library(Rtsne)
library(ggplot2)
library(RColorBrewer) #カラーパレット
library(dplyr)


# データのロード ---------------------------------------------------------------------

mat = fread(
  "PCAWG_matrix_6type.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)
labels = fread(
  "PCAWG_matrix_labels.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)

x1 = fread(
  "PCAWG_matrix_position.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)

gene = read.table(
  "TableS3_panorama_driver_mutations_ICGC_samples.public.tsv",
  sep = "",
  header = TRUE
)

label = as.vector(as.matrix(labels[, 1]))
type = as.vector(as.matrix(labels[, 2]))
barcode = as.vector(as.matrix(labels[, 3]))

# 距離行列 --------------------------------------------------------------------

x2 = matrix(rep(0, sum(nrow(mat)) * 6), nrow = nrow(mat))
for (i in 1:1933) {
  x2[i, ] = c(sum(mat[i, c(1:3127)]),
              sum(mat[i, c(3128:6254)]),
              sum(mat[i, c(6255:9381)]),
              sum(mat[i, c(9382:12508)]),
              sum(mat[i, c(12509:15635)]),
              sum(mat[i, c(15636:18762)]))
}

x1 = as.matrix(x1)
tic()
X2 <- x1 %*% t(x1)
d1 <- matrix(0, 1933, 1933)
for (i1 in 1:1933) {
  for (i2 in 1:1933) {
    d1[i1, i2] <- X2[i1, i1] - 2 * X2[i1, i2] + X2[i2, i2]
  }
}
X2 <- x2 %*% t(x2)
d2 <- matrix(0, 1933, 1933)
for (i1 in 1:1933) {
  for (i2 in 1:1933) {
    d2[i1, i2] <- X2[i1, i1] - 2 * X2[i1, i2] + X2[i2, i2]
  }
}
toc()


# geneの行列作成 ---------------------------------------------------------------
gene = as.matrix(gene)
x3 <- matrix(0, length(barcode), length(table(gene[, 7])))
colnames(x3) = sort(unique(gene[, 7]))
sorted = sort(unique(gene[, 7]))
j = 1
for (i in 1:length(barcode)) {
  sub = subset(gene, gene[, 1] %in% barcode[i])
  if (nrow(sub) != 0) {
    for (j in 1:nrow(sub)) {
      n = match(sub[j, 7], sorted)
      x3[i, n] = x3[i, n] + 1
    }
    
  }
}
X2 <- x3 %*% t(x3)
d3 <- matrix(0, 1933, 1933)
for (i1 in 1:1933) {
  for (i2 in 1:1933) {
    d3[i1, i2] <- X2[i1, i1] - 2 * X2[i1, i2] + X2[i2, i2]
  }
}

# tSNE --------------------------------------------------------------------
toPoint = function(factors) {
  mapping <- c (
    "Biliary-AdenoCA" = 1,
    "Bone-Cart" = 2,
    "Bone-Epith" = 3,
    "Bone-Osteosarc" = 4,
    "Breast-AdenoCa" = 5,
    "Breast-DCIS" = 6,
    "Breast-LobularCa" = 7,
    "CNS-Medullo" = 8,
    "CNS-PiloAstro" = 9,
    "Eso-AdenoCa" = 10,
    "Head-SCC" = 11,
    "Kidney-RCC" = 12,
    "Liver-HCC" = 13,
    "Lymph-BNHL" = 14,
    "Lymph-CLL" = 15,
    "Lymph-NOS" = 16,
    "Myeloid-AML" = 17,
    "Myeloid-MDS" = 18,
    "Myeloid-MPN" = 19,
    "Ovary-AdenoCA" = 20,
    "Panc-AdenoCA" = 21,
    "Panc-Endocrine" = 22,
    "Prost-AdenoCA" = 23,
    "Skin-Melanoma" = 24,
    "Stomach-AdenoCA" = 25
  )
  mapping[as.character(factors)]
}
type_num = as.integer(unlist(lapply(type, toPoint)))

# w1 = c(1, 0.75, 0.5, 0.25, 0)
# w2 = c(0, 0.25, 0.5, 0.75, 1)
w1 = c(1, 0, 0, 1/3)
w2 = c(0, 1, 0, 1/3)
w3 = c(0, 0, 1, 1/3)
for (i in 1:4) {
  D = w1[i] * d1 + w2[i] * d2 + w3[i] * d3
  tic()
  tsne = Rtsne(
    D,
    check_duplicates = FALSE,
    verbose = TRUE,
    initial_dims = nrow(D),
    is_distance = TRUE
  )
  toc()
  file = sprintf("~/Genome/tsne_2matrix_w1_%s_w2_%s_w3_%s.png",
                 w1[i],
                 w2[i],
                 w3[i])
  title = sprintf("tsne_2matrix_w1_%s_w2_%s_w3_%s.png", w1[i], w2[i], w3[i])
  png(file,
      width = 2000,
      height = 2000,)
  plot(tsne$Y, t = 'n', main = title)
  # legend("bottomleft", legend = sort(unique(type)), col = c(1:25),pch = c(1:25)%/%6+15)
  legend(
    "bottomleft",
    legend = sort(unique(type)),
    col = c(1:25),
    pch = c(1:25)
  )
  # text(tsne$Y, labels = as.character(as.factor(label)),col = type)
  # points(tsne$Y,col = type_num, pch = (type_num%/%6)+15)
  points(tsne$Y, col = type_num, pch = type_num)
  dev.off()
}
