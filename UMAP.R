library(data.table)
library(tictoc)
library(umap)
library(ggplot2)
library(RColorBrewer) #カラーパレット
library(dplyr) #カラーパレット


mat = fread(
  "matrix_6types.csv",
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

label = as.vector(as.matrix(labels[,1]))
type = as.vector(as.matrix(labels[,2]))


# UMAP --------------------------------------------------------------------

tic()
genome.umap = umap(mat)
toc()

#プロット関連
toPoint = function(factors) { 
  mapping <- c ("Biliary-AdenoCA"= 1,
                "Bone-Cart"= 2,
                "Bone-Epith"= 3,
                "Bone-Osteosarc"= 4,
                "Breast-AdenoCa"= 5,
                "Breast-DCIS"= 6,
                "Breast-LobularCa"= 7,
                "CNS-Medullo"= 8,
                "CNS-PiloAstro"= 9,
                "Eso-AdenoCa"= 10,
                "Head-SCC"= 11,
                "Kidney-RCC"= 12,
                "Liver-HCC"= 13,
                "Lymph-BNHL"= 14,
                "Lymph-CLL"= 15,
                "Lymph-NOS"= 16,
                "Myeloid-AML"= 17,
                "Myeloid-MDS"= 18,
                "Myeloid-MPN"= 19,
                "Ovary-AdenoCA"= 20,
                "Panc-AdenoCA"= 21,
                "Panc-Endocrine"= 22,
                "Prost-AdenoCA"= 23,
                "Skin-Melanoma"= 24,
                "Stomach-AdenoCA"= 25)
  mapping[as.character(factors)]
}
type_num = as.integer(unlist(lapply(type,toPoint)))
type_n = as.integer(type_num%%6+1)
k = unique(type)
s = c()
for (i in 1:25) {
  s = unlist(c(s,k[i]))
}
s=sort(s)
colPal3 <- colorRampPalette(brewer.pal(11, "Spectral"))
cp3 = colPal3(25)
#ここまで
g.layout=as.data.frame(genome.umap$layout)
g.type = as.data.frame(cbind(genome.umap$layout,type,as.factor(type_n)))
colnames(g.layout) = c("x","y")
colnames(g.type) = c("x","y","type","type_n")

png("~/Genome/umap_alldata_6type2.png",
    width = 2000,
    height = 2000,
)

print(
  ggplot(data=g.layout,aes(x=x,y=y))+
    geom_point(aes(color=factor(type),shape=type_n))+
    scale_shape_identity()+
    # geom_point(aes(color=type),shape=type_n)+
    ggtitle("UMAP for 25 cancers(6type)")
)

# print(
#   ggplot(g.type,aes(x,y))+
#     geom_point()+
#     scale_shape_identity()+
#     # geom_point(aes(color=type),shape=type_n)+
#     ggtitle("UMAP for 25 cancers(6type)")
# )
dev.off()  
