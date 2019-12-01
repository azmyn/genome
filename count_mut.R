library(data.table)
library(tictoc)


# データ読み込み -----------------------------------------------------------------


args <- commandArgs(trailingOnly = TRUE) #引数受け取り
tic()
# mc3 = fread(
#   "/Users/azumi/Genome/mc3/mc3.v0.2.8.PUBLIC.maf.gz",
#   stringsAsFactors = FALSE,
#   encoding = "UTF-8"
# )
# clinical = fread(
#   "/Users/azumi/Genome/laml_tcga/data_bcr_clinical_data_patient.txt",
#   stringsAsFactors = FALSE,
#   encoding = "UTF-8",
#   skip = 4
# )

mc3 = fread(args[1], stringsAsFactors = FALSE, encoding = "UTF-8")
clinical = fread(
  args[2],
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  skip = 4
)
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
# 変異数のカウント ----------------------------------------------------------------


sample_mut =  data.frame(matrix(rep(NA, 3), nrow = 1))[numeric(0), ]
colnames(sample_mut) = c("PATIENT_ID", "MUTATIONS", "AMOUNT_OF_SUBSTITUTIONS")
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

for (i in 1:length(data$Tumor_Sample_Barcode)) {
  newID = data$Tumor_Sample_Barcode[i]
  if (newID != ID) {
    j = j + 1
    ID = newID
  }
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
    # }else if(data$Variant_Type[i] == "INS"){
    #   m = paste0( "[",data$Reference_Allele[i],">",data$Tumor_Seq_Allele2[i] ,"]")
    # }else if(data$Variant_Type[i] == "DEL"){
    #   m = paste0( "[",data$Reference_Allele[i],">",data$Tumor_Seq_Allele2[i] ,"]")
    # }else if(data$Variant_Type[i] == "TNP"){
    #
    # }else if(data$Variant_Type[i] == "ONP"){
    #
    # }else if(data$Variant_Type[i] == "INS"||data$Variant_Type[i] == "DEL"||data$Variant_Type[i] == "ONP"||data$Variant_Type[i] == "TNP"||data$Variant_Type[i] == "DNP"){
  # } else if (data$Variant_Type[i] == "DNP") {
  #   m = paste0("[",
  #              data$Reference_Allele[i],
  #              ">",
  #              data$Tumor_Seq_Allele2[i] ,
  #              "]")
  #   
  } else{
    next
  }
  
  if (is.na(sample_mut[j, 1])) {
    sample_mut[j, 1] = as.character(substr(ID, 1, 12))
    sample_mut[j, 2] = m  
    sample_mut[j, 3] = 1
  } else{
    sample_mut[j, 2] = paste(sample_mut[j, 2], m)
    sample_mut[j, 3] = sample_mut[j, 3] + 1
  }
  
}
print(nrow(sample_mut))
print("個のサンプルごとに集計完了")
# 書き出し --------------------------------------------------------------------

file.name = sprintf("/Users/azumi/Dropbox/KU/shimolab_2019/genome/%s.mut_count.txt", substr(args[2], 21, 24))

fwrite(sample_mut, file = file.name, row.names =
         F) # 一度書き出し
