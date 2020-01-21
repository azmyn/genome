library(data.table)
library(tictoc)

# データ読み込み -----------------------------------------------------------------
tic()
data = fread(
  "/Users/azumi/Genome/f_c_p_screened.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)
toc()
# 変異数のカウント ----------------------------------------------------------------

sample_chrm_mut =  data.frame(matrix(rep(NA, 1), nrow = 1))[numeric(0), ]
sample_chrm_mut_type =  data.frame(matrix(rep(NA, 1), nrow = 1))[numeric(0), ]
sample_chrm_mut_pos =  data.frame(matrix(rep(NA, 1), nrow = 1))[numeric(0), ]
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

c = cbind(sample_chrm_mut, sample_chrm_mut_type, sample_chrm_mut_pos, )
colnames(c) = c("Mutations", "Mut_type", "Pos_1Mb")

toc() #だいたい10000秒ぐらい

d = as.data.table(cbind(data, as.data.table(c)))[, c(1, 4, 9, 11:16)]
fwrite(d, file = "mutation_small.csv", row.names = F) # 一度書き出し

# 行列作成(書き換えが必要) --------------------------------------------------------------------

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
TSB = unique(d$Tumor_Sample_Barcode)
mat = data.frame(matrix(rep(NA, sum(chrmlen_1Mb)*length(TSB)), ncol = sum(chrmlen_1Mb) ,nrow = length(TSB)))
label = unique(mutate(d, label = paste(Project_Code ,Tumor_Sample_Barcode,sep="_")) [,9])
TSB = unique(d$Tumor_Sample_Barcode)
type = c()
Donor_ID = c()
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
for (i in 1:length(TSB)) {
  sub = filter(d, Tumor_Sample_Barcode == TSB[i]) 
  type = c(type, sub$Project_Code[1])
  Donor_ID = c(Donor_ID, sub$Donor_ID[1])
  for (j in 1:nrow(sub)) {
    c = calc(sub$Chromosome[j], sub$Pos_1Mb[j])
    if (is.na(mat[i, c])) {
      mat[i, c] = 1
    } else{
      mat[i, c] = mat[i, c] + 1
    }
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
fwrite(as.data.table(cbind(label, type, TSB,Donor_ID)), file = "matrix_labels.csv", row.names = F)

toc() #10000秒ぐらい

# 変異別を書き直したぶん -----------------------------------------------------------------

mutation = c("C>A","C>G","C>T","T>A","T>C","T>G")

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
d = as.matrix(d)
for (i in 1:6) {
  mat_mut = matrix(0,length(barcode),sum(chrmlen_1Mb))
  tic(i)
  d_s = subset(d,d[,7] == mutation[i])
  for (j in 1:length(barcode)) {
    d_s_b = subset(d_s,d_s[,3]==barcode[j])
    for (k in 1:nrow(d_s_b)) {
      c = calc(d_s_b[k,1], d_s_b[k,8])
      if (is.na(mat_mut[j, c])) {
        mat_mut[j, c] = 1
      } else{
        mat_mut[j, c] = mat_mut[j, c] + 1
      }
    }
  }
  filename = sprintf("PCAWG_matrix_6type_part%s.csv",i)
  assign(paste("mat_mut_", i, sep=""), mat_mut)
  fwrite(mat_mut, filename, row.names = F)
  toc(log=TRUE)
}
mat_all = cbind(mat_mut_1,mat_mut_2,mat_mut_3,mat_mut_4,mat_mut_5,mat_mut_6)
fwrite(mat_all,"PCAWG_matrix_6type.csv",row.names = F)

x2 = matrix(rep(0, sum(nrow(mat)) * 6), nrow = nrow(mat))
for (i in 1:length(barcode)) {
  x2[i,] = c(sum(mat_all[i, c(1:3127)]),
             sum(mat_all[i, c(3128:6254)]),
             sum(mat_all[i, c(6255:9381)]),
             sum(mat_all[i, c(9382:12508)]),
             sum(mat_all[i, c(12509:15635)]),
             sum(mat_all[i, c(15636:18762)]))
}
fwrite(x2,"PCAWG_matrix_type.csv",row.names = F)
