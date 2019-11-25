library(data.table)
library(ggplot2)
library(ggbiplot)
library(tm)

# データ読み込み -----------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE) #引数受け取り 1つ目がmutation, 2つ目がclinical
mut_dat = fread(args[1], stringsAsFactors = FALSE, encoding = "UTF-8")
cli_dat = fread(
  args[2],
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  skip = 4
)

# 対応付け ---------------------------------------------------------------------

dat = cbind(mut_dat, data.frame(matrix(rep(NA, nrow(
  mut_dat
)), nrow = nrow(mut_dat))))


colnames(dat) = c("PATIENT_ID",
                  "MUTATIONS",
                  "AMOUNT_OF_MUTATIONS",
                  "SMOKING_HISTORY")

for (i in 1:length(mut_dat$PATIENT_ID)) {
  for (j in 1:length(cli_dat$PATIENT_ID)) {
    if (substr(mut_dat$PATIENT_ID[i], 1, 12) == cli_dat$PATIENT_ID[j]) {
      smoke = cli_dat$TOBACCO_SMOKING_HISTORY_INDICATOR[j]
      if (smoke == '[Not Available]') {
        dat$SMOKING_HISTORY[i] = NA
      } else if (smoke == 1) {
        dat$SMOKING_HISTORY[i] = "Lifelong Non-Smoker"
      } else if (smoke == 2) {
        dat$SMOKING_HISTORY[i] = "Current Smoker"
      } else if (smoke == 3) {
        dat$SMOKING_HISTORY[i] = "Current Reformed Smoker for > 15 yrs"
      } else if (smoke == 4) {
        dat$SMOKING_HISTORY[i] = "Current Refomed Smoker for < or = 15 yrs"
      } else if (smoke == 5) {
        dat$SMOKING_HISTORY[i] = "Current Reformed Smoker Duration Not Specified"
        #二値の場合
        # } else if (smoke == 2 ||
        #            smoke == 3 || smoke == 4 || smoke == 5) {
        #   mut_dat[i] = "SMOKER"
      } else {
        dat$SMOKING_HISTORY[i] = NA
        
      }
    }
  }
}

dat = na.omit(dat[, 1:4])

# 頻度行列 --------------------------------------------------------------------

corpus = Corpus(VectorSource(dat$MUTATIONS))
tdm = TermDocumentMatrix(corpus)
D = as.matrix(tdm)
rownames(D) = toupper(rownames(D))

# PCA ---------------------------------------------------------------------

dpca = prcomp(t(D))

file.name = sprintf("%s_biplot.png", substr(args[1], 1, 4))
png(file.name, width = 2000, height = 2000)
print(
  ggbiplot(
    dpca,
    obs.scale = 1,
    var.scale = 1,
    groups = dat$SMOKING_HISTORY,
    ellipse = TRUE,
    circle = TRUE
  )
)
dev.off()
