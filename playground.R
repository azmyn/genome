mat = fread(
  "/Users/azumi/Genome/matrix.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)

fwrite(as.list(table(labels$type)),"temp.csv",row.names = F)
labels= fread("matrix_labels.csv", stringsAsFactors = FALSE,
              encoding = "UTF-8",
              sep = ","
)
table(labels$type)



df_f3 = subset(df_f, labels$type %in% c("Breast-AdenoCa","Breast-DCIS","Breast-LobularCa"))

data = subset(data, data$FILTER %in% "PASS")


X <- dataGen(m=1, n=100, p=10, eps=0.2, bLength=4)$data[[1]]

resR <- robpca(X, k=2)
diagPlot(resR)