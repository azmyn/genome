library(data.table)
library(tictoc)


# データ読み込み -----------------------------------------------------------------


# args <- commandArgs(trailingOnly = TRUE) #引数受け取り
tic()
mc3 = fread(
  "/Users/azumi/Genome/mc3/mc3.v0.2.8.PUBLIC.maf.gz",
  stringsAsFactors = FALSE,
  encoding = "UTF-8"
)
clinical = fread(
  "/Users/azumi/Genome/hnsc_tcga/data_bcr_clinical_data_patient.txt",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  skip = 4
)

# mc3 = fread(args[1], stringsAsFactors = FALSE, encoding = "UTF-8")
# clinical = fread(
#   args[2],
#   stringsAsFactors = FALSE,
#   encoding = "UTF-8",
#   skip = 4
# )
toc()

print("データ読み込み完了")
# mc3から該当する種類のデータの抽出 ------------------------------------------------------

getid <- function(x) {
  return(substr(x, 1, 12))
}
barcode = numeric(0)

barcode = apply(as.matrix(mc3$Tumor_Sample_Barcode), 1, getid)
data = subset(mc3, barcode %in% clinical$PATIENT_ID)
data = subset(data, data$Variant_Type %in% "SNP") #SNPのみ抽出
data = subset(data, data$FILTER %in% "PASS")

print(nrow(data))
print("行の該当データのみ抽出")


# 変異のカウント(mutationのvariationを考慮しない )---------------------------------------


sample_chrm_mut = data[, c(16, 5, 6)]
sample_chrm_mut = cbind(sample_chrm_mut,  data.frame(matrix(
  rep(1,  nrow(data)), nrow = nrow(data), ncol = 1
)))
colnames(sample_chrm_mut) = c(
  "PATIENT_ID",
  "CHROMOSOME",
  "START_POSITION",
  # "MUTATIONS",
  "AMOUNT_OF_SUBSTITUTIONS"
)

# 1Mbごとに切り分け --------------------------------------------------------------

sorted = sample_chrm_mut[with(sample_chrm_mut,
                              order(PATIENT_ID, CHROMOSOME, START_POSITION)),]

mut_1Mb = data.frame(matrix(rep(NA, 4), nrow = 1))[numeric(0), ]
colnames(mut_1Mb) = c("PATIENT_ID",
                      "CHROMOSOME",
                      "POSITION",
                      # "MUTATIONS",
                      "AMOUNT_OF_SUBSTITUTIONS")

j = 1
chrm = sorted[1, 2]
tic()
for (i in 1:nrow(sorted[, 1])) {
  ID = sorted[i, 1]
  newchrm = sorted[i, 2]
  pos = sorted[i, 3] %/% 1000000
  if (newchrm != chrm) {
    j = j + 1
    chrm = newchrm
    mut_1Mb[j, 1] = ID
    mut_1Mb[j, 2] = newchrm
    mut_1Mb[j, 3] = pos
    # mut_1Mb[j, 4] = sorted[i, 4]
    mut_1Mb[j, 4] = 1
  } else if (i == 1) {
    mut_1Mb[j, 1] = ID
    mut_1Mb[j, 2] = newchrm
    mut_1Mb[j, 3] = pos
    # mut_1Mb[j, 4] = sorted[i, 4]
    mut_1Mb[j, 4] = 1
  } else if (pos != mut_1Mb[j, 3]) {
    j = j + 1
    mut_1Mb[j, 1] = ID
    mut_1Mb[j, 2] = newchrm
    mut_1Mb[j, 3] = pos
    # mut_1Mb[j, 4] = sorted[i, 4]
    mut_1Mb[j, 4] = 1
  } else if(chrm == newchrm && ID != sorted[i-1, 1]){
    j = j + 1
    mut_1Mb[j, 1] = ID
    mut_1Mb[j, 2] = newchrm
    mut_1Mb[j, 3] = pos
  } else{
    # mut_1Mb[j, 4] = paste(mut_1Mb[j, 4], sorted[i, 4])
    mut_1Mb[j, 4] = mut_1Mb[j, 4] + 1
  }
}
toc()  #ここ2時間ぐらいかかる

print(nrow(mut_1Mb))
print("個のサンプルごとに集計完了")
#  ここ!!!!!!!!
# file.name = sprintf("~/Genome/%s.mut_1Mb_count.png", args[3])
file.name = "~/Genome/hnsc.mut_1Mb_count.txt"
fwrite(mut_1Mb, file = file.name, row.names =
         F) # 一度書き出し
