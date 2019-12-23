library(data.table)
library(tictoc)

# データ読み込み -----------------------------------------------------------------


# args <- commandArgs(trailingOnly = TRUE) #引数受け取り
tic()
# data = fread("/Users/azumi/Genome/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz", stringsAsFactors = FALSE, encoding = "UTF-8")
data = fread("/Users/azumi/Dropbox/KU/shimolab_2019/genome/test.csv", stringsAsFactors = FALSE, encoding = "UTF-8", sep=",")

toc()

print("データ読み込み完了")
# data2 =data
dat = data[1:10000,]
# 変異数のカウント ----------------------------------------------------------------

sample_chrm_mut =  data.frame(matrix(rep(NA, 1), nrow = 1))[numeric(0), ]
# colnames(sample_chrm_mut) = c(
#   "PATIENT_ID",
#   "CHROMOSOME",
#   "START_POSITION",
#   "MUTATIONS",
#   "AMOUNT_OF_SUBSTITUTIONS"
# )
# i = 1
# j = 1
# ID = data$Tumor_Sample_Barcode[1]

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
  # newID = data$Tumor_Sample_Barcode[i]
  # if (newID != ID) {
  #   j = j + 1
  #   ID = newID
  # }
  # if (data$Variant_Type[i] == "SNP" &&
  #     data$FILTER[i] == "PASS" && data$NCALLERS[i] >= 3) {
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

  # if (is.na(sample_chrm_mut[j, 1])) {
  # sample_chrm_mut[j, 1] = data$Donor_ID[i]
  # sample_chrm_mut[j, 2] = data$Chromosome[i]
  # sample_chrm_mut[j, 3] = data$Start_position[i]
  sample_chrm_mut[i] = toupper(m)
  # sample_chrm_mut[j, 5] = data$Project_Code[i]
  # }
  # } else{
  #   sample_chrm_mut[j, 2] = data$Chromosome[i]
  #   sample_chrm_mut[j, 3] = data$Start_Position[i]
  #   sample_chrm_mut[j, 4] = paste(sample_chrm_mut[j, 4], m)
  #   # sample_chrm_mut[j, 5] = sample_chrm_mut[j, 5] + 1
  # }
  #
  # j = j + 1
}
toc()
print(nrow(sample_chrm_mut))
print("個のサンプルと染色体ごとに集計完了")

fwrite(as.list(sample_chrm_mut), file = "mutation.csv", row.names = F) # 一度書き出し
         

# あ -----------------------------------------------------------------------
tic()

sorted = sample_chrm_mut[with(sample_chrm_mut,
order(PATIENT_ID, CHROMOSOME, START_POSITION)),]

mut_1Mb = data.frame(matrix(rep(NA, 5), nrow = 1))[numeric(0), ]
colnames(mut_1Mb) = c("PATIENT_ID",
                      "CHROMOSOME",
                      "POSITION",
                      "MUTATIONS",
                      "AMOUNT_OF_SUBSTITUTIONS")

j = 1
chrm = sorted[1, 2]

for (i in 1:nrow(sorted[, 1])) {
  ID = sorted[i, 1]
  newchrm = sorted[i, 2]
  pos = sorted[i, 3] %/% 1000000
  mut = sorted[i, 4]
  if (newchrm != chrm) {
    j = j + 1
    chrm = newchrm
    mut_1Mb[j, 1] = ID
    mut_1Mb[j, 2] = newchrm
    mut_1Mb[j, 3] = pos
    mut_1Mb[j, 4] = mut
    mut_1Mb[j, 5] = 1
  } else if (i == 1) {
    mut_1Mb[j, 1] = ID
    mut_1Mb[j, 2] = newchrm
    mut_1Mb[j, 3] = pos
    mut_1Mb[j, 4] = mut
    mut_1Mb[j, 5] = 1
  } else if (pos != mut_1Mb[j, 3]) {
    j = j + 1
    mut_1Mb[j, 1] = ID
    mut_1Mb[j, 2] = newchrm
    mut_1Mb[j, 3] = pos
    mut_1Mb[j, 4] = mut
    mut_1Mb[j, 5] = 1
  } else{
    mut_1Mb[j, 4] = paste(mut_1Mb[j, 4], sorted[i, 4])
    mut_1Mb[j, 5] = mut_1Mb[j, 5] + 1
  }
}
toc()  #ここ2時間ぐらいかかる

print(nrow(mut_1Mb))
print("個のサンプルごとに集計完了")

# file.name = sprintf("~/Genome/%s.mut_1Mb_count.txt", args[3])
file.name = "all.mut_1Mb_count.txt"
fwrite(mut_1Mb, file = file.name, row.names =
         F) # 一度書き出し
