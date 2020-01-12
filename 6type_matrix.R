library(data.table)
library(tictoc)

d1 = fread(
  "mutation_small.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)

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

# 1爪 ----------------------------------------------------------------------


d = subset(d1, d1$Mut_type %in% "C>G")
mat = data.frame(matrix(rep(0, sum(chrmlen_1Mb)), nrow = 1))[numeric(0),]
ID =  d$Donor_ID[1]
j = 1
label = sprintf("%s_%s", d$Project_Code[1], d$Donor_ID[1])
type = d$Project_Code[1]

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

for (i in 1:nrow(d)) {
  newID = d$Donor_ID[i]
  if (newID != ID) {
    j = j + 1
    ID = newID
    label = c(label, sprintf("%s_%s", d$Project_Code[i], d$Donor_ID[i]))
    type = c(type, d$Project_Code[i])
  }
  # if(d$Mut_type[i] == "C>A"){
  c = calc(d$Chromosome[i], d$Pos_1Mb[i])
  if (is.na(mat[j, c])) {
    mat[j, c] = 1
  } else{
    mat[j, c] = mat[j, c] + 1
  }
  # }
  
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
    n = sprintf("%s-%s-C>G", c, as.character(j))
    cname = c(cname, as.vector(n))
  }
}
colnames(mat) = cname
# rownames(mat) = label
fwrite(mat, file = "matrix2_6.csv", row.names = F) # 一度書き出し

toc()
mat2 = mat


# もう一個 --------------------------------------------------------------------

d = subset(d1, d1$Mut_type %in% "C>T")
mat = data.frame(matrix(rep(0, sum(chrmlen_1Mb)), nrow = 1))[numeric(0),]
ID =  d$Donor_ID[1]
j = 1
label = sprintf("%s_%s", d$Project_Code[1], d$Donor_ID[1])
type = d$Project_Code[1]

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

for (i in 1:nrow(d)) {
  newID = d$Donor_ID[i]
  if (newID != ID) {
    j = j + 1
    ID = newID
    label = c(label, sprintf("%s_%s", d$Project_Code[i], d$Donor_ID[i]))
    type = c(type, d$Project_Code[i])
  }
  # if(d$Mut_type[i] == "C>A"){
  c = calc(d$Chromosome[i], d$Pos_1Mb[i])
  if (is.na(mat[j, c])) {
    mat[j, c] = 1
  } else{
    mat[j, c] = mat[j, c] + 1
  }
  # }
  
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
    n = sprintf("%s-%s-C>T", c, as.character(j))
    cname = c(cname, as.vector(n))
  }
}
colnames(mat) = cname
# rownames(mat) = label
fwrite(mat, file = "matrix3_6.csv", row.names = F) # 一度書き出し
toc()
mat3=mat

# もう一個 --------------------------------------------------------------------

d = subset(d1, d1$Mut_type %in% "T>A")
mat = data.frame(matrix(rep(0, sum(chrmlen_1Mb)), nrow = 1))[numeric(0),]
ID =  d$Donor_ID[1]
j = 1
label = sprintf("%s_%s", d$Project_Code[1], d$Donor_ID[1])
type = d$Project_Code[1]

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

for (i in 1:nrow(d)) {
  newID = d$Donor_ID[i]
  if (newID != ID) {
    j = j + 1
    ID = newID
    label = c(label, sprintf("%s_%s", d$Project_Code[i], d$Donor_ID[i]))
    type = c(type, d$Project_Code[i])
  }
  # if(d$Mut_type[i] == "C>A"){
  c = calc(d$Chromosome[i], d$Pos_1Mb[i])
  if (is.na(mat[j, c])) {
    mat[j, c] = 1
  } else{
    mat[j, c] = mat[j, c] + 1
  }
  # }
  
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
    n = sprintf("%s-%s-T>A", c, as.character(j))
    cname = c(cname, as.vector(n))
  }
}
colnames(mat) = cname
# rownames(mat) = label
fwrite(mat, file = "matrix4_6.csv", row.names = F) # 一度書き出し
toc()
mat4=mat

# もう一個 --------------------------------------------------------------------

d = subset(d1, d1$Mut_type %in% "T>C")
mat = data.frame(matrix(rep(0, sum(chrmlen_1Mb)), nrow = 1))[numeric(0),]
ID =  d$Donor_ID[1]
j = 1
label = sprintf("%s_%s", d$Project_Code[1], d$Donor_ID[1])
type = d$Project_Code[1]

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

for (i in 1:nrow(d)) {
  newID = d$Donor_ID[i]
  if (newID != ID) {
    j = j + 1
    ID = newID
    label = c(label, sprintf("%s_%s", d$Project_Code[i], d$Donor_ID[i]))
    type = c(type, d$Project_Code[i])
  }
  # if(d$Mut_type[i] == "C>A"){
  c = calc(d$Chromosome[i], d$Pos_1Mb[i])
  if (is.na(mat[j, c])) {
    mat[j, c] = 1
  } else{
    mat[j, c] = mat[j, c] + 1
  }
  # }
  
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
    n = sprintf("%s-%s-T>C", c, as.character(j))
    cname = c(cname, as.vector(n))
  }
}
colnames(mat) = cname
# rownames(mat) = label
fwrite(mat, file = "matrix5_6.csv", row.names = F) # 一度書き出し
toc()
mat5=mat

# もう一個 --------------------------------------------------------------------

d = subset(d1, d1$Mut_type %in% "T>G")
mat = data.frame(matrix(rep(0, sum(chrmlen_1Mb)), nrow = 1))[numeric(0),]
ID =  d$Donor_ID[1]
j = 1
label = sprintf("%s_%s", d$Project_Code[1], d$Donor_ID[1])
type = d$Project_Code[1]

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

for (i in 1:nrow(d)) {
  newID = d$Donor_ID[i]
  if (newID != ID) {
    j = j + 1
    ID = newID
    label = c(label, sprintf("%s_%s", d$Project_Code[i], d$Donor_ID[i]))
    type = c(type, d$Project_Code[i])
  }
  # if(d$Mut_type[i] == "C>A"){
  c = calc(d$Chromosome[i], d$Pos_1Mb[i])
  if (is.na(mat[j, c])) {
    mat[j, c] = 1
  } else{
    mat[j, c] = mat[j, c] + 1
  }
  # }
  
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
    n = sprintf("%s-%s-T>G", c, as.character(j))
    cname = c(cname, as.vector(n))
  }
}
colnames(mat) = cname
# rownames(mat) = label
fwrite(mat, file = "matrix6_6.csv", row.names = F) # 一度書き出し
toc()
mat6=mat

mat1 = fread(
  "matrix1_6.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)

merge(mat1, mat2,by) -> newmat 
newmat = cbind(mat1, mat2,mat3, mat4,mat5,mat6)

fwrite(newmat,file="matrix_6types.csv",row.names = F)
table(newmat[1,])
