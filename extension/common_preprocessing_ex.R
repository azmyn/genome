library(data.table)
library(tictoc)
library(dplyr)
library(tidyverse)
# データ読み込み -----------------------------------------------------------------
tic()
data = fread(
  "/Users/azumi/Genome/f_c_p_screened_ex.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)
toc() #80秒ぐらい
# 変異数のカウント ----------------------------------------------------------------

sample_chrm_mut = c()
sample_chrm_mut_type = c()
sample_chrm_mut_pos = c()

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
DNP_R = c("GT", "GG", "AG", "GA", "CA", "AA")
DNP_COMP = c("AT", "CG", "GC", "TA")
DNP_AT = c("CC", "GA", "CA")
DNP_CG = c("GT", "TC", "TT")
DNP_GC = c("AA", "AG", "CA")
DNP_TA = c("CT", "GG", "GT")

tic()
for (i in 1:nrow(data)) {
  if (data$Variant_Type[i] == "SNP") {
    ref = as.character(data$Reference_Allele[i])
    if (ref == "T" || ref == "C") {
      fb_minus =  str_sub(data$ref_context[i], 10, 10)
      fb_plus =  str_sub(data$ref_context[i], 12, 12)
      mut = as.character(data$Tumor_Seq_Allele2[i])
    } else{
      fb_minus =  comp(str_sub(data$ref_context[i], 10, 10))
      fb_plus =  comp(str_sub(data$ref_context[i], 12, 12))
      ref = comp(ref)
      mut = comp(as.character(data$Tumor_Seq_Allele2[i]))
    }
    m = paste0(fb_minus, "[", ref, ">", mut, "]", fb_plus)
    m_type = paste0(ref, ">", mut)
    sample_chrm_mut_type[i] = toupper(m_type)
  } else if (data$Variant_Type[i] == "DNP") {
    if (data$Reference_Allele[i] %in% DNP_COMP) {
      if (data$Tumor_Seq_Allele2[i] %in% DNP_COMP) {
        ref_a = data$Reference_Allele[i]
        tum_a = data$Tumor_Seq_Allele2[i]
      } else{
        if(data$Reference_Allele[i] == "AT" && data$Tumor_Seq_Allele2[i]%in% DNP_AT){
          ref_a = data$Reference_Allele[i]
          tum_a = data$Tumor_Seq_Allele2[i]
        }else if (data$Reference_Allele[i] == "CG" && data$Tumor_Seq_Allele2[i]%in% DNP_CG){
          ref_a = data$Reference_Allele[i]
          tum_a = data$Tumor_Seq_Allele2[i]
        }else if (data$Reference_Allele[i] == "GC" && data$Tumor_Seq_Allele2[i]%in% DNP_GC){
          ref_a = data$Reference_Allele[i]
          tum_a = data$Tumor_Seq_Allele2[i]
        }else if (data$Reference_Allele[i] == "TA" && data$Tumor_Seq_Allele2[i]%in% DNP_TA){
          ref_a = data$Reference_Allele[i]
          tum_a = data$Tumor_Seq_Allele2[i]
        }else{
          ref_a = data$Reference_Allele[i]
          tum_a = paste0(comp(str_sub(data$Tumor_Seq_Allele2[i], 2, 2)), comp(str_sub(data$Tumor_Seq_Allele2[i], 1, 1)))
        }
      }
    } else if (data$Reference_Allele[i] %in% DNP_R) {
      ref_a = paste0(comp(str_sub(data$Reference_Allele[i], 2, 2)), comp(str_sub(data$Reference_Allele[i], 1, 1)))
      tum_a = paste0(comp(str_sub(data$Tumor_Seq_Allele2[i], 2, 2)), comp(str_sub(data$Tumor_Seq_Allele2[i], 1, 1)))
    } else{
      ref_a = data$Reference_Allele[i]
      tum_a = data$Tumor_Seq_Allele2[i]
    }
    m = paste0("[", ref_a, ">", tum_a , "]")
    sample_chrm_mut_type[i] = toupper(m) #変異タイプごとに6種類
  } else if (data$Variant_Type[i] == "INS" ||
             data$Variant_Type[i] == "DEL") {
    ref_a = data$Reference_Allele[i]
    tum_a = data$Tumor_Seq_Allele2[i]
    m = paste0("[", ref_a, ">", tum_a , "]")
    sample_chrm_mut_type[i] = toupper(m) #変異タイプごとに6種類
  } else{
    print(i)
    print("error!")
  }
  sample_chrm_mut[i] = toupper(m) #前後の文脈を含めて96種類に分類
  # sample_chrm_mut_type[i] = toupper(m_type) #変異タイプごとに6種類
  sample_chrm_mut_pos[i] = data$Start_position[i] %/% 1000000 # 1Mで割った商
}

c = cbind(sample_chrm_mut, sample_chrm_mut_type, sample_chrm_mut_pos)
colnames(c) = c("Mutations", "Mut_type", "Pos_1Mb")

toc() #だいたい6500秒ぐらい

d = as.data.table(cbind(data, as.data.table(c)))[, c(1, 5, 9, 11:15)]
fwrite(d, file = "mutation_small_ex.csv", row.names = F) # 一度書き出し



# やっぱりSNPとDNPだけにする --------------------------------------------------------

d = filter(d, Variant_Type == "DNP"|Variant_Type == "SNP")

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
TSB = unique(d$Tumor_Sample_Barcode)

mat = data.frame(matrix(
  rep(NA, sum(chrmlen_1Mb) * 1),
  ncol = sum(chrmlen_1Mb) ,
  nrow = 1
))

label = unique(mutate(d, label = paste(Project_Code , Tumor_Sample_Barcode, sep =
                                         "_")) [, 9])
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
fwrite(mat, file = "PCAWG_matrix_ex.csv", row.names = F) # 一度書き出し
fwrite(as.data.table(cbind(label, type, TSB, Donor_ID)), file = "PCAWG_matrix_labels_ex.csv", row.names = F)

toc() #8000秒ぐらい

# 文脈を含めたSNPとDNPの174パターン ------------------------------------------------------------

mutation = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
nu = c("A", "C", "G", "T")
mutation96 = c()
for (i in 1:6) {
  for (j in 1:4) {
    for (k in 1:4) {
      mutation96 = c(mutation96, paste0(nu[j], "[", mutation[i], "]", nu[k]))
    }
  }
}
dnp =  filter(d, Variant_Type == "DNP")
mutation174 = c (mutation96, sort(unique(dnp$Mutations)))

mat_174 =  data.frame(matrix(rep(NA, 1), nrow = 1, ncol = 174))[numeric(0),]

d = as.data.frame(d)
tic()
for (i in 1:length(TSB)) {
  sub = filter(d, d$Tumor_Sample_Barcode == TSB[i])
  for (j in 1:nrow(sub)) {
    c = match(sub$Mutations[j], mutation174)
    if (is.na(mat_174[i, c])) {
      mat_174[i, c] = 1
    } else{
      mat_174[i, c] = mat_174[i, c] + 1
    }
  }
}
mat_174[is.na(mat_174)] <- 0
toc() #3000秒ぐらい
fwrite(mat_174, "PCAWG_matrix_context_174type_ex.csv", row.names = F)
