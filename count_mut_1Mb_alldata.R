library(data.table)
library(tictoc)

# データ読み込み -----------------------------------------------------------------


# args <- commandArgs(trailingOnly = TRUE) #引数受け取り
tic()
# data = fread("/Users/azumi/Genome/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz", stringsAsFactors = FALSE, encoding = "UTF-8")
data = fread(
  "/Users/azumi/Genome/f_c_p_screened.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)
toc()

print("データ読み込み完了")

dat = data[1:20000, ]
# 変異数のカウント ----------------------------------------------------------------

sample_chrm_mut =  data.frame(matrix(rep(NA, 1), nrow = 1))[numeric(0),]
sample_chrm_mut_type =  data.frame(matrix(rep(NA, 1), nrow = 1))[numeric(0),]
sample_chrm_mut_pos =  data.frame(matrix(rep(NA, 1), nrow = 1))[numeric(0),]

comp <- function(x) {
  #変異を揃える
  x = toupper(x)
  switch(
    x,
    "A" = return("T"),
    "T" = return("A"),
    "G" = return("C"),
    "C" = return("G"),
    print("?")
  )
}
tic()
for (i in 1:nrow(data)) {
  ref = as.character(data$Reference_Allele[i])
  if (ref == "T" || ref == "C") {
    fb_minus =  substr(data$ref_context[i], 10, 10)
    fb_plus =  substr(data$ref_context[i], 12, 12)
    mut = as.character(data$Tumor_Seq_Allele2[i])
  } else{
    fb_minus =  comp(substr(data$ref_context[i], 10, 10))
    fb_plus =  comp(substr(data$ref_context[i], 12, 12))
    ref = comp(ref)
    mut = comp(as.character(data$Tumor_Seq_Allele2[i]))
  }
  m = paste0(fb_minus, "[", ref, ">", mut, "]", fb_plus)
  m_type = paste0(ref, ">", mut)
  sample_chrm_mut[i] = toupper(m) #前後の文脈を含めて96種類に分類
  sample_chrm_mut_type[i] = toupper(m_type) #変異タイプごとに6種類
  sample_chrm_mut_pos[i] = data$Start_position[i] %/% 1000000 # 1Mで割った商
}
print(nrow(sample_chrm_mut))
print("個のサンプルと染色体ごとに集計完了")

c = cbind(sample_chrm_mut, sample_chrm_mut_type, sample_chrm_mut_pos)
colnames(c) = c("Mutations", "Mut_type", "Pos_1Mb")

toc() #だいたい10000秒ぐらい

d = as.data.table(cbind(data, as.data.table(c)))[, c(1, 4, 11:15)]
fwrite(d, file = "mutation_small.csv", row.names = F) # 一度書き出し

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
  # 57227415 Y染色体の数がおかしい
  58227415
)
chrmlen_1Mb = ceiling(chrmlen / 1000000 + 1)

mat = data.frame(matrix(rep(0, sum(chrmlen_1Mb)), nrow = 1))[numeric(0),]
ID =  d$Donor_ID[1]
j = 1
label = sprintf("%s_%s", d$Project_Code[1], d$Donor_ID[1])
type = d$Project_Code[1]

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

for (i in 1:nrow(d)) {
  newID = d$Donor_ID[i]
  if (newID != ID) {
    j = j + 1
    ID = newID
    label = c(label, sprintf("%s_%s", d$Project_Code[i], d$Donor_ID[i]))
    type = c(type, d$Project_Code[i])
  }
  c = calc(d$Chromosome[i], d$Pos_1Mb[i])
  if (is.na(mat[j, c])) {
    mat[j, c] = 1
  } else{
    mat[j, c] = mat[j, c] + 1
  }
}

mat[is.na(mat)] <- 0

cname = NULL
for (i in 1:length(chrmlen_1Mb)) {
  if (i == 23) {
    c = 'X'
  } else if (i == 24) {
    c = 'Y'
  } else{
    c = i
  }
  for (j in 1:chrmlen_1Mb[i]) {
    n = sprintf("%s-%s", c, as.character(j))
    cname = c(cname, as.vector(n))
  }
}
colnames(mat) = cname
# rownames(mat) = label
fwrite(mat, file = "matrix.csv", row.names = F) # 一度書き出し
fwrite(as.data.table(cbind(label, type)), file = "matrix_labels.csv", row.names = F) # 一度書き出し

# 変異別行列作成 -----------------------------------------------------------------

mat_mut = data.frame(matrix(rep(0, sum(chrmlen_1Mb) * 6), nrow = 1))[numeric(0),]
ID =  d$Donor_ID[1]
j = 1
label = sprintf("%s_%s", d$Project_Code[1], d$Donor_ID[1])
type = d$Project_Code[1]

calc = function(x, y, z) {
  if (z == "C>A") {
    mut = 1
  } else if (z == "C>G") {
    mut = 2
  } else if (z == "C>T") {
    mut = 3
  } else if (z == "T>A") {
    mut = 4
  } else if (z == "T>C") {
    mut = 5
  } else if (z == "T>G") {
    mut = 6
  } else{
    stop("error")
  }
  if (x == 'X') {
    x = 23
  }
  if (x == 'Y') {
    x = 24
  }
  if (x == 1) {
    return((mut - 1) * sum(chrmlen_1Mb) + y + 1)
  } else{
    x = as.numeric(x)
    return ((mut - 1) * sum(chrmlen_1Mb) + sum(chrmlen_1Mb[1:(x - 1)]) + y + 1)
  }
}
tic()

for (i in 1:nrow(d)) {
  newID = d$Donor_ID[i]
  if (newID != ID) {
    j = j + 1
    ID = newID
    label = c(label, sprintf("%s_%s", d$Project_Code[i], d$Donor_ID[i]))
    type = c(type, d$Project_Code[i])
  }
  c = calc(d$Chromosome[i], d$Pos_1Mb[i], d$Mut_type[i])
  if (is.na(mat_mut[j, c])) {
    mat_mut[j, c] = 1
  } else{
    mat_mut[j, c] = mat_mut[j, c] + 1
  }
}

mat_mut[is.na(mat_mut)] <- 0
toc()