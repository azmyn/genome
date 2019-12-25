library(data.table)
library(tictoc)
library(ggbiplot)
library(rpca)

mat = fread(
  "matrix.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)
labels = fread(
  "matrix_labels.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)

label = labels[,1]
type = labels[,2]


# 行列(確率ベクトルバージョン) ---------------------------------------------------------

sum = apply(mat, 1, sum)
prob <- function(x) {
  return(x / sum(x))
}
mat2 = as.matrix (mat)
mat_prob = matrix(rep(0, sum(chrmlen_1Mb) * nrow(mat2)), nrow = nrow(mat2))
for (i in 1:nrow(mat)) {
  for (j in 1:ncol(mat)) {
    mat_prob[i, j] = mat2[i, j] / sum[i]
  }
}
colnames(mat_prob) = cname

# PCA ---------------------------------------------------------------------

# df = as.data.frame(mat)
df = as.data.frame(mat_prob)
df_f <- df[, apply(df, 2, var, na.rm = TRUE) != 0]


df_f2 = df_f[c(-1288,-1763),]
dpca = prcomp(df_f2)
tic()
dpca = robpca(df_f,k=2)
toc()
# biplot(x=dpca)

# file.name = sprintf("/Users/azumi/Genome/%s_and_%s.1Mb_PCA.png", args[1], args[2])
# title = sprintf("%s(red) and %s(green)", toupper(args[1]), toupper(args[2]))

file.name = "alldata_prob_screened2.png"
title = "All 25 cancer"

png(file.name,
    width = 2000,
    height = 2000)
print(
  ggbiplot(
    dpca,
    obs.scale = 1,
    var.scale = 1,
    groups = type[c(-1288,-1763)],
    ellipse = TRUE,
    circle = TRUE,
    color = type[c(-1288,-1763)],
  )
  + ggtitle(title)
)
dev.off()


# PCA_small ---------------------------------------------------------------

# df = as.data.frame(mat_prob[c(-712,-1288),])
df = as.data.frame(mat_prob[c(-1863,-1879),])
df3 = subset(df, labels$type %in% c("Breast-AdenoCa","Breast-DCIS","Breast-LobularCa"))[c(-109,-121),]
# df3 = subset(df, labels$type %in% c("Lymph-BNHL","Lymph-CLL"))[c(-85,-144),]
df_f3 <- df3[, apply(df3, 2, var, na.rm = TRUE) != 0]
label_small = subset(labels$type, labels$type %in% c("Breast-AdenoCa","Breast-DCIS","Breast-LobularCa"))[c(-109,-121)]
# label_small = subset(labels$type,labels$type %in% c("Lymph-BNHL","Lymph-CLL"))[c(-85,-144)]

dpca = prcomp(df_f3)
biplot(x=dpca)

file.name = "B_cancer_data_prob_.png"
title = "Breast cancer(3types)"

png(file.name,
    width = 2000,
    height = 2000)
print(
  ggbiplot(
    dpca,
    obs.scale = 1,
    var.scale = 1,
    groups = label_small,
    ellipse = TRUE,
    circle = TRUE,
    color = label_small
  )
  + ggtitle(title)
)
dev.off()