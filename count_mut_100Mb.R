library(data.table)
library(tictoc)


# データ読み込み -----------------------------------------------------------------


args <- commandArgs(trailingOnly = TRUE) #引数受け取り
tic()
mc3 = fread(
  "/Users/azumi/Genome/mc3/mc3.v0.2.8.PUBLIC.maf.gz",
  stringsAsFactors = FALSE,
  encoding = "UTF-8"
)
clinical = fread(
  "/Users/azumi/Genome/lusc_tcga/data_bcr_clinical_data_patient.txt",
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

print(nrow(data))
print("行の該当データのみ抽出")

data = subset(data, data$FILTER %in% "PASS")



# 変異数のカウント ----------------------------------------------------------------


sample_chrm_mut =  data.frame(matrix(rep(NA, 5), nrow = 1))[numeric(0), ]
colnames(sample_chrm_mut) = c(
  "PATIENT_ID",
  "CHROMOSOME",
  "START_POSITION",
  "MUTATIONS",
  "AMOUNT_OF_SUBSTITUTIONS"
)
i = 1
j = 1
ID = data$Tumor_Sample_Barcode[1]

comp <- function(x) {
  #変異を揃える
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
  if (data$Variant_Type[i] == "SNP" &&
      data$FILTER[i] == "PASS" && data$NCALLERS[i] >= 3) {
    # if(data$Variant_Type[i] == "SNP") {
    ref = as.character(data$Reference_Allele[i])
    if (ref == "T" || ref == "C") {
      fb_minus =  substr(data$CONTEXT[i], 5, 5)
      fb_plus =  substr(data$CONTEXT[i], 7, 7)
      mut = as.character(data$Tumor_Seq_Allele2[i])
    } else{
      fb_minus =  comp(substr(data$CONTEXT[i], 5, 5))
      fb_plus =  comp(substr(data$CONTEXT[i], 7, 7))
      ref = comp(ref)
      mut = comp(as.character(data$Tumor_Seq_Allele2[i]))
    }
    m = paste0(fb_minus, "[", ref, ">", mut, "]", fb_plus)
    
  } else{
    next
  }
  
  # if (is.na(sample_chrm_mut[j, 1])) {
  sample_chrm_mut[j, 1] = as.character(substr(data$Tumor_Sample_Barcode[i], 1, 12))
  sample_chrm_mut[j, 2] = data$Chromosome[i]
  sample_chrm_mut[j, 3] = data$Start_Position[i]
  sample_chrm_mut[j, 4] = m
  sample_chrm_mut[j, 5] = 1
  # } else{
  #   sample_chrm_mut[j, 2] = data$Chromosome[i]
  #   sample_chrm_mut[j, 3] = data$Start_Position[i]
  #   sample_chrm_mut[j, 4] = paste(sample_chrm_mut[j, 4], m)
  #   # sample_chrm_mut[j, 5] = sample_chrm_mut[j, 5] + 1
  # }
  #
  j = j + 1
}
toc()
print(nrow(sample_chrm_mut))
print("個のサンプルと染色体ごとに集計完了")


# 変異のカウント(mutationのvariationを考慮しない )---------------------------------------
nrow(mc3)

sample_chrm_mut = data[, c(16, 5, 6)]
sample_chrm_mut = cbind(sample_chrm_mut,  data.frame(matrix(
  rep(1,  nrow(data)), nrow = nrow(data), ncol = 2
)))
colnames(sample_chrm_mut) = c(
  "PATIENT_ID",
  "CHROMOSOME",
  "START_POSITION",
  "MUTATIONS",
  "AMOUNT_OF_SUBSTITUTIONS"
)
# 書き出し --------------------------------------------------------------------

file.name = sprintf(
  "/Users/azumi/Dropbox/KU/shimolab_2019/genome/%s.chrm_mut_count.txt",
  substr(args[2], 21, 24)
) #がんの略称が4文字じゃないと問題を起こすので注意

fwrite(sample_chrm_mut, file = file.name, row.names =
         F) # 一度書き出し

sample_chrm_mut = fread("/Users/azumi/Dropbox/KU/shimolab_2019/genome/luad.chrm_mut_count.txt")
# 1Mbごとに切り分け --------------------------------------------------------------

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
tic()
for (i in 1:nrow(sorted[, 1])) {
  newchrm = sorted[i, 2]
  pos = sorted[i, 3] %/% 1000000
  if (newchrm != chrm) {
    j = j + 1
    chrm = newchrm
    mut_1Mb[j, 1] = sorted[i, 1]
    mut_1Mb[j, 2] = sorted[i, 2]
    mut_1Mb[j, 3] = pos
    mut_1Mb[j, 4] = sorted[i, 4]
    mut_1Mb[j, 5] = 1
  } else if (i == 1) {
    mut_1Mb[j, 1] = sorted[i, 1]
    mut_1Mb[j, 2] = sorted[i, 2]
    mut_1Mb[j, 3] = pos
    mut_1Mb[j, 4] = sorted[i, 4]
    mut_1Mb[j, 5] = 1
  } else if (pos != mut_1Mb[j, 3]) {
    j = j + 1
    mut_1Mb[j, 1] = sorted[i, 1]
    mut_1Mb[j, 2] = sorted[i, 2]
    mut_1Mb[j, 3] = pos
    mut_1Mb[j, 4] = sorted[i, 4]
    mut_1Mb[j, 5] = 1
  } else{
    mut_1Mb[j, 4] = paste(mut_1Mb[j, 4], sorted[i, 4])
    mut_1Mb[j, 5] = mut_1Mb[j, 5] + 1
  }
}
toc()

print(nrow(mut_1Mb))
print("個のサンプルごとに集計完了")

file.name = "lusc.mut_1Mb_count.txt"
fwrite(mut_1Mb, file = file.name, row.names =
         F) # 一度書き出し
