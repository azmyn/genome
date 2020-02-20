library(data.table)
library(tictoc)
# ロード ---------------------------------------------------------------------
setwd("~/Genome/PCAWG")

pos_vec = fread(
  "position_lda_vector.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)
type_vec = fread(
  "type96_lda_vector.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)
gene_vec = fread(
  "gene_lda_vector.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)
labels = fread(
  "PCAWG_matrix_labels.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)


label = as.vector(as.matrix(labels[, 1]))
type = as.vector(as.matrix(labels[, 2]))
barcode = as.vector(as.matrix(labels[, 3]))
donor = as.vector(as.matrix(labels[, 4]))

# -------------------------------------------------------------------------

posRF = as.data.frame(pos_vec[,1:50])
typeRF = as.data.frame(type_vec[,1:50])
geneRF = as.data.frame(gene_vec[,1:30])
attach(posRF)


no.forests=25 # for the final version,you would want to increase this number to say 50 or 100
no.trees=3000 

tic()
pos_distRF = RFdist(posRF, mtry1=3, no.trees, no.forests, addcl1=T,addcl2=F,imp=T, oob.prox1=T)
toc()
tic()
type_distRF = RFdist(typeRF, mtry1=3, no.trees, no.forests, addcl1=T,addcl2=F,imp=T, oob.prox1=T)
toc()
tic()
gene_distRF = RFdist(geneRF, mtry1=3, no.trees, no.forests, addcl1=T,addcl2=F,imp=T, oob.prox1=T)
toc()

tic()
no.clusters = 2
pos_labelRF = pamNew(pos_distRF$cl1, no.clusters)
toc()
tic()
pos_labelEuclid = pamNew(dist(as.data.frame(posRF)), no.clusters)
toc()
tic()
type_labelRF = pamNew(type_distRF$cl1, no.clusters)
toc()
tic()
type_labelEuclid = pamNew(dist(as.data.frame(typeRF)), no.clusters)
toc()
tic()
gene_labelRF = pamNew(gene_distRF$cl1, no.clusters)
toc()
tic()
gene_labelEuclid = pamNew(dist(as.data.frame(geneRF)), no.clusters)
toc()


table(pos_labelRF)

no.clusters = 2
labelRF = pamNew(distRF$cl1, no.clusters)
labelEuclid = pamNew(dist(as.data.frame(posRF)), no.clusters)
fisher.test(table(labelRF, labelEuclid))

labelNew = ifelse(labelRF==1&labelEuclid==1, 1, 
                  ifelse(labelRF==1&labelEuclid==2, 2,
                         ifelse(labelRF==2&labelEuclid==1, 3, 4)))
fit1 = survfit(Surv(time, event)~labelNew, data=dat1, conf.type="log-log")
mylegend=c("RF cluster 1, Euclid cluster 1", "RF cluster 1, Euclid cluster 2",  
           "RF cluster 2, Euclid cluster 1","RF cluster 2, Euclid cluster 2")
plot(fit1, conf.int=F,col= unique(labelNew), lty=1:4, xlab="Time to death ",ylab="Survival",legend.text=mylegend, lwd=1,mark.time=TRUE) 


# 以下テストデータ ----------------------------------------------------------------


setwd("~/Dropbox/KU/shimolab_2019/genome/appendix")
source("FunctionsRFclustering.txt")
## read in the data set
## This is the data set we used in the technical report Shi and Horvath (2005)
## as the motivational example
## We will show how to generate the plots of Figure 1 in that manuscript
dat1 = read.table("testData.csv", sep=",", header=T, row.names=1)
## This is the input file for RF clustering algorithm
datRF = dat1[,1:8]
attach(datRF)
## Here is the histogram of tumor marker #1 as shown in Figure 1a
hist(datRF$Marker1, xlim=c(0,100), ylim=c(0,300), xlab="Score in %", main="Marker 1")
## Calculating RF distance between samples based on the 8 marker measurements
## This will take quite long time depends how many tree and repetitions you choose
## We suggest to use relatively large number of forests with large number of trees
no.forests=25 # for the final version,you would want to increase this number to say 50 or 100
no.trees=3000 # this could also be increased to say 4000
# Since we are mainly interested in the Addcl1 RF dissimilarity we set addcl1=T,addcl2=F
# imp=T specificies that we are also interested in the importance measures.
distRF = RFdist(datRF, mtry1=3, no.trees, no.forests, addcl1=T,addcl2=F,imp=T, oob.prox1=T) 

## PAM clustering based on the Addcl1 RF dissimilarity
no.clusters = 2
labelRF = pamNew(distRF$cl1, no.clusters)
## PAM clustering based on Euclidean distance
labelEuclid = pamNew(dist(datRF), no.clusters)
## Due to the randomness of RF procedure, the exact distance measure will vary a bit
## Therefore, we also include our RF clustering result in our data
##If you want to see our result, you may need to add the following statement
labelRF = dat1$labelRF
## Check the agreement between RF cluster and Euclidean distance cluster
fisher.test(table(labelRF, labelEuclid)) ## Fisher’s exact p value
# Selected Output
# Fishers Exact Test for Count Data
# data: table(labelRF, labelEuclid)
# p-value = 1.216e-15
## Define a new clustering label based on labelRF and labelEuclid
## labelNew=1, if labelRF=1 and labelEuclid=1
## labelNew=2, if labelRF=1 and labelEuclid=2
## labelNew=3, if labelRF=2 and labelEuclid=1
## labelNew=4, if labelRF=2 and labelEuclid=2
labelNew = ifelse(labelRF==1&labelEuclid==1, 1,
 ifelse(labelRF==1&labelEuclid==2, 2,
 ifelse(labelRF==2&labelEuclid==1, 3, 4)))
## check survival difference as in Figure 1b
## variables "time" and "event" in dat1 are survival time and cencering indicator, respectively
## NOTE: the RF clusters are more meaningful with respect to survival time.
fit1 = survfit(Surv(time, event)~labelNew, data=dat1, conf.type="log-log")
mylegend=c("RF cluster 1, Euclid cluster 1", "RF cluster 1, Euclid cluster 2",
 "RF cluster 2, Euclid cluster 1","RF cluster 2, Euclid cluster 2")
plot(fit1, conf.int=F,col= unique(labelNew), lty=1:4, xlab="Time to death
",ylab="Survival",legend.text=mylegend, lwd=1,mark.time=TRUE) 


#THE END

