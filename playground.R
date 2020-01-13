library(data.table)
library(tictoc)
d = fread(
  "mutation_small.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)
labels =  fread(
  "PCAWG_matrix_labels.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)

label = as.vector(as.matrix(labels[, 1]))
type = as.vector(as.matrix(labels[, 2]))
barcode = as.vector(as.matrix(labels[, 3]))

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
  tic()
  d_s = subset(d,d$Mut_type %in% mutation[i])
  for (j in 1:length(barcode)) {
    d_s_b = subset(d_s,d_s$Tumor_Sample_Barcode %in% barcode[j])
    for (k in 1:nrow(d_s_b)) {
      c = calc(d_s_b$Chromosome[k], d_s_b$Pos_1Mb[k])
      if (is.na(mat_mut[j, c])) {
        mat_mut[j, c] = 1
      } else{
        mat_mut[j, c] = mat_mut[j, c] + 1
      }
    }
    filename = sprintf("PCAWG_matrix_6type_part%s.csv",i)
    assign(paste("mat_mut_", i, sep=""), mat_mut)
    fwrite(mat_mut, filename, row.names = F) # 一度書き出し
  }
  toc(log=TRUE)
}


