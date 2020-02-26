library(data.table)
library(tictoc)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(gtable)

# データのロード ---------------------------------------------------------------------
setwd("~/Genome/PCAWG_ex")

labels = fread(
  "PCAWG_matrix_labels_ex.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)

x1 = fread(
  "PCAWG_matrix_position_ex.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)
# x2 = fread(
#   "PCAWG_matrix_type.csv",
#   stringsAsFactors = FALSE,
#   encoding = "UTF-8",
#   sep = ","
# )
x2 = fread(
  "PCAWG_matrix_context_174type.csv",
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

# 距離行列 --------------------------------------------------------------------

x1 = as.matrix(x1)
x2 = as.matrix(x2)
X1 = x1 %*% t(x1)
X2 = x2 %*% t(x2)
X3 = x3 %*% t(x3)

d1 = matrix(0, length(barcode), length(barcode))
d2 = matrix(0, length(barcode), length(barcode))
d3 = matrix(0, length(barcode), length(barcode))

for (i1 in 1:length(barcode)) {
  for (i2 in 1:length(barcode)) {
    d1[i1, i2] = sqrt(X1[i1, i1] - 2 * X1[i1, i2] + X1[i2, i2])
    d2[i1, i2] = sqrt(X2[i1, i1] - 2 * X2[i1, i2] + X2[i2, i2])
    d3[i1, i2] = sqrt(X3[i1, i1] - 2 * X3[i1, i2] + X3[i2, i2])
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
toShape = function(factors) {
  mapping <- c (
    "Biliary-AdenoCA" = 0,
    "Bone-Cart" = 1,
    "Bone-Epith" = 1,
    "Bone-Osteosarc" = 1,
    "Breast-AdenoCa" = 2,
    "Breast-DCIS" = 2,
    "Breast-LobularCa" = 2,
    "CNS-Medullo" = 3,
    "CNS-PiloAstro" = 3,
    "Eso-AdenoCa" = 4,
    "Head-SCC" = 5,
    "Kidney-RCC" = 6,
    "Liver-HCC" = 7,
    "Lymph-BNHL" = 8,
    "Lymph-CLL" = 8,
    "Lymph-NOS" = 8,
    "Myeloid-AML" = 9,
    "Myeloid-MDS" = 9,
    "Myeloid-MPN" = 9,
    "Ovary-AdenoCA" = 10,
    "Panc-AdenoCA" = 11,
    "Panc-Endocrine" = 11,
    "Prost-AdenoCA" = 12,
    "Skin-Melanoma" = 13,
    "Stomach-AdenoCA" = 14
  )
  mapping[as.character(factors)]
}
shape_num = as.integer(unlist(lapply(type, toShape)))
tic()
w1 = c(1, 0, 0, 0.5, 0.5, 0, 1 / 3)
w2 = c(0, 1, 0, 0.5, 0, 0.5, 1 / 3)
w3 = c(0, 0, 1, 0, 0.5, 0.5, 1 / 3)
perp = seq(5, 50, by = 5)
# perp = 30
for (i in 1:length(w1)) {
  for (j in 1:length(perp)) {
    D = w1[i] * d1 + w2[i] * d2 + w3[i] * d3
    tsne = Rtsne(
      D,
      check_duplicates = FALSE,
      verbose = TRUE,
      initial_dims = nrow(D),
      is_distance = TRUE,
      perplexity = perp[j]
    )
    file = sprintf(
      "~/Genome/PCAWG_ex/3matrix/tsne_3matrix_perp_%s_w1_%s_w2_%s_w3_%s.png",
      perp[j],
      signif(w1[i], digits = 3),
      signif(w2[i], digits = 3),
      signif(w3[i], digits = 3)
    )
    title = sprintf(
      "SNP/DNP 3matrix perp=%s w1=%s w2=%s w3=%s",
      perp[j],
      signif(w1[i], digits = 3),
      signif(w2[i], digits = 3),
      signif(w3[i], digits = 3)
    )
    df <- data.frame(matrix(rep(NA, 3), nrow = 1950))[numeric(0), ]
    df = as.data.frame(cbind(as.factor(type), tsne$Y[, 1] , tsne$Y[, 2]))
    df[, 1] = as.factor(type)
    df[, 4] = as.factor(shape_num)
    colnames(df) <- c("gene_type", "tSNE_1", "tSNE_2", "pch")
    g <-
      ggplot(df,
             aes(
               x = df$tSNE_1,
               y = df$tSNE_2,
               color = df$gene_type,
               shape = df$pch
             )) +
      geom_point() +
      scale_color_manual(
        name = "Cancer detail types",
        labels = sort(unique(type)),
        values = c(
          "#F5D174",
          "#F3C0AB",
          "#E07987",
          "#D45D87",
          "#E7A5C9",
          "#CF8CBB",
          "#BB9CD2",
          "#AEC1E3",
          "#5FAFD7",
          "#75D4C9",
          "#8DDA81",
          "#CADF77",
          "#F5D174",
          "#F3C0AB",
          "#E07987",
          "#D45D87",
          "#E7A5C9",
          "#CF8CBB",
          "#BB9CD2",
          "#AEC1E3",
          "#5FAFD7",
          "#75D4C9",
          "#8DDA81",
          "#CADF77",
          "#DA5019"
        )
      ) +
      scale_shape_manual(
        name = "Cancer types",
        labels = c(
          "Biliary",
          "Bone",
          "Breast",
          "CNS",
          "Eso",
          "Head",
          "Kidney",
          "Liver",
          "Lymph",
          "Myeloid",
          "Ovary",
          "Panc",
          "Prost",
          "Skin",
          "Stomach"
        ),
        values = c(0:14)
      ) +
      labs(x = "tSNE_1", y = "tSNE_2") +
      ggtitle(title)
    ggsave(
      file = file,
      plot = g,
      dpi = 500,
      width = 16,
      height = 9
    )
  }
}

toc()