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

