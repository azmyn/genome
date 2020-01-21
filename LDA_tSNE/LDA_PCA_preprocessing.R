library(data.table)
library(tictoc)
library(tm)
library(NLP)
library(lda)
library(dplyr)

# ロード ---------------------------------------------------------------------
setwd("~/Genome/PCAWG")

d = fread(
  "mutation_small.csv",
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
label = as.vector(as.matrix(labels[, 1]))
type = as.vector(as.matrix(labels[, 2]))
barcode = as.vector(as.matrix(labels[, 3]))

# geneの行列作成 ---------------------------------------------------------------

x3 <- matrix(0, length(barcode), length(table(gene[, 7])))
colnames(x3) = sort(unique(gene[, 7]))
sorted = sort(unique(gene[, 7]))
j = 1
for (i in 1:length(barcode)) {
  sub = filter(gene, sample_id == barcode[i])
  if (nrow(sub) != 0) {
    for (j in 1:nrow(sub)) {
      n = match(sub[j, 7], sorted)
      x3[i, n] = x3[i, n] + 1
    }
  }
}

# 文章の形にする ----------------------------------------------------------------------
d2 = as.data.frame(d)
d2 = mutate(d2, NewCol = paste(!!!rlang::syms(c(
  "Chromosome", "Pos_1Mb"
)), sep = "_"))

sent = matrix(NA, length(barcode), 3)
tic()
for (i in 1:length(barcode)) {
  sub = filter(d2, Tumor_Sample_Barcode ==  barcode[i])
  for (j in 1:nrow(sub)) {
    if (is.na(sent[i, 1])) {
      sent[i, 1] = sub[j, 9]
    } else{
      sent[i, 1] = paste(sent[i, 1], sub[j, 9])
    }
    if (is.na(sent[i, 2])) {
      sent[i, 2] = as.character(sub[j, 7])
    } else{
      sent[i, 2] = paste(sent[i, 2], as.character(sub[j, 7]))
    }
  }
  for (k in 1:ncol(x3)) {
    if (x3[i, k] != 0) {
      if (is.na(sent[i, 3])) {
        sent[i, 3] = as.character(k)
      } else{
        sent[i, 3] = paste(sent[i, 3], as.character(k))
      }
    }
  }
  print(i)
}
fwrite(sent, "sent.csv", row.names = F)
toc()#70000秒ぐらい


# LDA ---------------------------------------------------------------------

sent = as.matrix(sent)
position_lex = lexicalize(sent[, 1], lower = FALSE)
type_lex = lexicalize(sent[, 2], lower = FALSE)
gene_lex = lexicalize(sent[, 3], lower = FALSE)
tic()
topic = 50
topic_gene = 30
pos_result <-
  lda.collapsed.gibbs.sampler(
    position_lex$documents,
    topic,
    position_lex$vocab,
    100,
    0.1,
    0.1,
    compute.log.likelihood = TRUE
  )
type_result <-
  lda.collapsed.gibbs.sampler(type_lex$documents,
                              topic,
                              type_lex$vocab,
                              100,
                              0.1,
                              0.1,
                              compute.log.likelihood = TRUE)
gene_result <-
  lda.collapsed.gibbs.sampler(
    gene_lex$documents,
    topic_gene,
    gene_lex$vocab,
    100,
    0.1,
    0.1,
    compute.log.likelihood = TRUE
  )
toc()

pos_vec = matrix(NA, length(label), topic)
type_vec = matrix(NA, length(label), topic)
gene_vec = matrix(NA, length(label), topic_gene)

for (i in 1:length(label)) {
  for (j in 1:topic) {
    pos_vec[i, j] = pos_result$document_sums[j, i] / sum(pos_result$document_sums[, i])
    type_vec[i, j] = type_result$document_sums[j, i] / sum(type_result$document_sums[, i])
  }
  for (j in 1:topic_gene) {
    if (sum(gene_result$document_sums[, i]) == 0) {
      gene_vec[i, j] = 0
    } else{
      gene_vec[i, j] = gene_result$document_sums[j, i] / sum(gene_result$document_sums[, i])
    }
  }
}
fwrite(pos_vec, "position_lda_vector.csv", row.names = F)
fwrite(type_vec, "type_lda_vector.csv", row.names = F)
fwrite(gene_vec, "gene_lda_vector.csv", row.names = F)
