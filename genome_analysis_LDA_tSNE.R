library(data.table)
library(tictoc)
library(tm)
library(NLP)
library(lda)
library(dplyr)
library(Rtsne)
library(ggplot2)
# library(RColorBrewer)#カラーパレット
library(tidytext)
library(stringr)
library(openxlsx)

# ロード ---------------------------------------------------------------------
setwd("~/Genome/PCAWG")

pos_vec = fread(
  "position_lda_vector.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)
type_vec = fread(
  "type_lda_vector.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)
gene_vec = fread(
  "gene_lda_vector.csv",
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

gene = read.table(
  "TableS3_panorama_driver_mutations_ICGC_samples.public.tsv",
  sep = "",
  header = TRUE
)

clinical_data_base = read.xlsx("pcawg_donor_clinical_August2016_v9.xlsx")

label = as.vector(as.matrix(labels[, 1]))
type = as.vector(as.matrix(labels[, 2]))
barcode = as.vector(as.matrix(labels[, 3]))
donor = as.vector(as.matrix(labels[, 4]))

# 距離行列 --------------------------------------------------------------------

d_pos = matrix(0, length(label), length(label))
d_type = matrix(0, length(label), length(label))
d_gene = matrix(0, length(label), length(label))
p = as.matrix(pos_vec) %*% t(as.matrix(pos_vec))
t = as.matrix(type_vec) %*% t(as.matrix(type_vec))
g = as.matrix(gene_vec) %*% t(as.matrix(gene_vec))

for (i in 1:length(label)) {
  for (j in 1:length(label)) {
    d_pos[i, j] = sqrt(p[i, i] - 2 * p[i, j] + p[j, j])
    d_type[i, j] = sqrt(t[i, i] - 2 * t[i, j] + t[j, j])
    d_gene[i, j] = sqrt(g[i, i] - 2 * g[i, j] + g[j, j])
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

w1 = c(1, 0, 0, 0.5, 0.5, 0, 1 / 3)
w2 = c(0, 1, 0, 0.5, 0, 0.5, 1 / 3)
w3 = c(0, 0, 1, 0, 0.5, 0.5, 1 / 3)
perp = seq(5, 50, by = 5)
# perp = 30
for (i in 1:length(w1)) {
  for (j in 1:length(perp)) {
    D = w1[i] * d_pos + w2[i] * d_type + w3[i] * d_gene
    tsne = Rtsne(
      D,
      check_duplicates = FALSE,
      verbose = TRUE,
      initial_dims = nrow(D),
      is_distance = TRUE,
      perplexity = perp[j],
    )
    file = sprintf(
      "~/Genome/LDA_tSNE/perplexity/lda_tsne_perp_%s_w1_%s_w2_%s_w3_%s.png",
      perp[j],
      signif(w1[i], digits = 3),
      signif(w2[i], digits = 3),
      signif(w3[i], digits = 3)
    )
    title = sprintf("LDA tSNE perplexity=%s w1=%s w2=%s w3=%s",
                    perp[j],
                    signif(w1[i], digits = 3),
                    signif(w2[i], digits = 3),
                    signif(w3[i], digits = 3)
                    )
    df <- data.frame(matrix(rep(NA, 3), nrow = 1950))[numeric(0),]
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


#clinical_smoking -------------------------------------------------------------------------
clinical_data = subset(clinical_data_base,
                       clinical_data_base$icgc_donor_id %in% donor)

smoking = filter(
  clinical_data,
  is.na(tobacco_smoking_history_indicator) == FALSE &
    tobacco_smoking_history_indicator != "Smoking history not documented"
)
compare = match(donor, smoking$icgc_donor_id)
comp = c()
class = c()
for (i in 1:1950) {
  if (is.na(compare[i]) == FALSE) {
    comp = c(comp, i)
    class = c (class, smoking$tobacco_smoking_history_indicator[compare[i]])
  }
}

ClasstoPoint = function(factors) {
  mapping <- c (
    "Lifelong non-smoker (<100 cigarettes smoked in lifetime)" = 5,
    "Current reformed smoker, duration not specified" = 4,
    "Current reformed smoker for > 15 years" = 3,
    "Current reformed smoker for <= 15 years"  = 2,
    "Current smoker (includes daily smokers non-daily/occasional smokers)" = 1
  )
  mapping[as.character(factors)]
}
class_num = as.integer(unlist(lapply(class, ClasstoPoint)))

w1 = c(1, 0, 0, 0.5, 0.5, 0, 1 / 3)
w2 = c(0, 1, 0, 0.5, 0, 0.5, 1 / 3)
w3 = c(0, 0, 1, 0, 0.5, 0.5, 1 / 3)

for (i in 1:length(w1)) {
  for (j in 1:length(perp)) {
    D = w1[i] * d_pos + w2[i] * d_type + w3[i] * d_gene
    tsne = Rtsne(
      D,
      check_duplicates = FALSE,
      verbose = TRUE,
      initial_dims = nrow(D),
      is_distance = TRUE,
    )
    file = sprintf("~/Genome/LDA_tSNE/smoking/smoking_perp_%s_w1_%s_w2_%s_w3_%s.png",
                   # perp[j],
                   signif(w1[i], digits = 3),
                   signif(w2[i], digits = 3),
                   signif(w3[i], digits = 3)
                   )
    title = sprintf("Smoking history w1=%s w2=%s w3=%s",
                    # perp[j],
                    signif(w1[i], digits = 3),
                    signif(w2[i], digits = 3),
                    signif(w3[i], digits = 3)
                    )
    df = data.frame(matrix(rep(NA, 3), nrow = 1950))[numeric(0),]
    df = as.data.frame(cbind(as.factor(class), tsne$Y[comp, 1] , tsne$Y[comp, 2]))
    df[, 1] = as.factor(class)
    df[, 4] = as.factor(class_num)
    colnames(df) = c("gene_type", "tSNE_1", "tSNE_2", "pch")
    g <-
      ggplot(df, aes(
        x = df$tSNE_1,
        y = df$tSNE_2,
        color = df$gene_type
      )) +
      geom_point() +
      scale_color_manual(
        name = "Smoking hisotry",
        labels = sort(unique(class)),
        values = c("#F3C0AB", "#D45D87", "#CF8CBB", "#5FAFD7", "#75D4C9")
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
#clinical_alcohol -------------------------------------------------------------------------
clinical_data = subset(clinical_data_base,
                       clinical_data_base$icgc_donor_id %in% donor)

alcohol = filter(clinical_data, alcohol_history == "yes" |
                   alcohol_history == "no")
compare = match(donor, alcohol$icgc_donor_id)
comp = c()
class = c()
for (i in 1:1950) {
  if (is.na(compare[i]) == FALSE) {
    comp = c(comp, i)
    class = c (class, alcohol$alcohol_history[compare[i]])
  }
}

ClasstoPoint = function(factors) {
  mapping <- c ("yes" = 1,
                "no" = 2)
  mapping[as.character(factors)]
}
class_num = as.integer(unlist(lapply(class, ClasstoPoint)))

w1 = c(1, 0, 0, 0.5, 0.5, 0, 1 / 3)
w2 = c(0, 1, 0, 0.5, 0, 0.5, 1 / 3)
w3 = c(0, 0, 1, 0, 0.5, 0.5, 1 / 3)
perp = seq(5, 50, by = 5)
perp = 30
for (i in 1:length(w1)) {
  for (j in 1:length(perp)) {
    D = w1[i] * d_pos + w2[i] * d_type + w3[i] * d_gene
    tsne = Rtsne(
      D,
      check_duplicates = FALSE,
      verbose = TRUE,
      initial_dims = nrow(D),
      is_distance = TRUE,
      perplexity = perp[j],
    )
    file = sprintf("~/Genome/LDA_tSNE/alcohol/alcohol_perp_%s_w1_%s_w2_%s_w3_%s.png",
                   # perp[j],
                   perp,
                   signif(w1[i], digits = 3),
                   signif(w2[i], digits = 3),
                   signif(w3[i], digits = 3)
    )
    title = sprintf("Alcohol history w1=%s w2=%s w3=%s",
                    # perp[j],
                    perp,
                    signif(w1[i], digits = 3),
                    signif(w2[i], digits = 3),
                    signif(w3[i], digits = 3)
    )
    df <- data.frame(matrix(rep(NA, 3), nrow = 1950))[numeric(0),]
    df = as.data.frame(cbind(as.factor(class), tsne$Y[comp, 1] , tsne$Y[comp, 2]))
    df[, 1] = as.factor(class)
    df[, 4] = as.factor(class_num)
    colnames(df) <- c("gene_type", "tSNE_1", "tSNE_2", "pch")
    g <-
      ggplot(df, aes(
        x = df$tSNE_1,
        y = df$tSNE_2,
        color = df$gene_type
      )) +
      geom_point() +
      scale_color_manual(
        name = "Alcohol history",
        labels = sort(unique(class)),
        values = c("#D45D87", "#5FAFD7")
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

# # subtype用にがん種を限定(データがなくてできなかった分, plotがggplot2になってないので注意) ------------------------------------------------------------------
# subtype_brca = subset(type, str_detect(type, "^Breast"))
# type_num = as.integer(unlist(lapply(subtype_brca, toPoint)))
# d_pos = d_type = d_gene= matrix(0, length(subtype_brca), length(subtype_brca))
#
# pos_brca = subset(pos_vec, str_detect(type, "^Breast"))
# type_brca = subset(type_vec, str_detect(type, "^Breast"))
# gene_brca = subset(gene_vec, str_detect(type, "^Breast"))
# p = as.matrix(pos_brca) %*% t(as.matrix(pos_brca))
# t = as.matrix(type_brca) %*% t(as.matrix(type_brca))
# g = as.matrix(gene_brca) %*% t(as.matrix(gene_brca))
#
# for (i in 1:length(subtype_brca)) {
#   for (j in 1:length(subtype_brca)) {
#     d_pos[i, j] = p[i, i] - 2 * p[i, j] + p[j, j]
#     d_type[i, j] = t[i, i] - 2 * t[i, j] + t[j, j]
#     d_gene[i, j] = g[i, i] - 2 * g[i, j] + g[j, j]
#   }
# }
#
# w1 = c(1, 0, 0, 0.5, 0.5, 0, 1 / 3)
# w2 = c(0, 1, 0, 0.5, 0, 0.5, 1 / 3)
# w3 = c(0, 0, 1, 0, 0.5, 0.5, 1 / 3)
# perp = seq(5, 50, by = 5)
# perp = 30
# for (i in 1:length(w1)) {
#   # for (j in 1:length(perp)) {
#   D = w1[i] * d_pos + w2[i] * d_type + w3[i] * d_gene
#   tsne = Rtsne(
#     D,
#     check_duplicates = FALSE,
#     verbose = TRUE,
#     initial_dims = nrow(D),
#     is_distance = TRUE,
#     # perplexity = perp[j],
#   )
#   file = sprintf("~/Genome/LDA_tSNE/breast_perp_%s_w1_%s_w2_%s_w3_%s.png",
#                  # perp[j],
#                  perp,
#                  w1[i],
#                  w2[i],
#                  w3[i])
#   title = sprintf("breast_perplexity_%s_w1_%s_w2_%s_w3_%s.png",
#                   # perp[j],
#                   perp,
#                   w1[i],
#                   w2[i],
#                   w3[i])
#
#   png(file,
#       width = 2000,
#       height = 2000,)
#   plot(tsne$Y, t = 'n', main = title)
#   legend(
#     "bottomleft",
#     legend = sort(unique(subtype_brca)),
#     col = c(1:3),
#     pch = c(1:3)
#   )
#   points(tsne$Y, col = type_num-4, pch = type_num-4)
#   dev.off()
#   # }
# }