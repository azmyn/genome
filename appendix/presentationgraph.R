graycolor = rep("#a9a9a9",25)
circle = rep(1,15)
oc = c(
  "#F5D174",
  "#F3C0AB",
  "#E07987",
  "#D45D87",
  "#E7A5C9",
  "#CF8CBB",
  "#BB9CD2",
  "#AEC1E3",
  "#5FAFD7",
  "#75D4C9",
  "#8DDA81",
  "#CADF77",
  "#F5D174",
  "#F3C0AB",
  "#E07987",
  "#D45D87",
  "#E7A5C9",
  "#CF8CBB",
  "#BB9CD2",
  "#AEC1E3",
  "#5FAFD7",
  "#75D4C9",
  "#8DDA81",
  "#CADF77",
  "#DA5019"
)
opch =  c(0:14)


# 肝臓がん胆道がん ----------------------------------------------------------------
# feature = "biliary_and_liver"
# fn = c(1,13)
# fgn = c(1,8)
# graycolor[1] = oc[25]

# 乳がん卵巣がん前立腺がん ------------------------------------------------------------
# feature = "breast_and_ovary_prost"
# fn = c(5,6,7,20,23)
# fgn = c(3,11,13)


# CNS ---------------------------------------------------------------------
# feature = "CNS"
# fn = c(8,9)
# fgn = c(4)


# 食道と胃がん ------------------------------------------------------------------
# feature = "eso_stomach"
# fn = c(10,25)
# fgn = c(5,15)


# リンパ 白血 ------------------------------------------------------------------
# feature = "lymph_myeloid"
# fn = c(14,15,16,17,18,19)
# fgn = c(9,10)

# 膵臓 ------------------------------------------------------------------
# feature = "pancreatic"
# fn = c(21,22)
# fgn = c(12)

# 骨 ------------------------------------------------------------------
feature = "bone"
fn = c(2,3,4)
fgn = c(2)


# 固定 ----------------------------------------------------------------------


graycolor[fn] = oc[fn]
circle[fgn] = opch[fgn]
circle[fgn] = opch[11]


# 3matrix -----------------------------------------------------------------


file = sprintf(
  "~/Genome/presentation/tsne_3matrix_%s_presen.png",
  feature
)


title = sprintf(
  "3matrix featuring : %s",
  feature
)
df <- data.frame(matrix(rep(NA, 3), nrow = 1950))[numeric(0), ]
df = as.data.frame(cbind(as.factor(type), tsne$Y[, 1] , tsne$Y[, 2]))
df[, 1] = as.factor(type)
df[, 4] = as.factor(shape_num)
colnames(df) <- c("gene_type", "tSNE_1", "tSNE_2", "pch")
g <-
  ggplot(df,
         aes(
           x = df$tSNE_1,
           y = df$tSNE_2,
           color = df$gene_type,
           shape = df$pch
         )) +
  geom_point() +
  scale_color_manual(
    name = "Cancer detail types",
    labels = sort(unique(type)),
    values = graycolor
  ) +
  scale_shape_manual(
    name = "Cancer types",
    labels = c(
      "Biliary",
      "Bone",
      "Breast",
      "CNS",
      "Eso",
      "Head",
      "Kidney",
      "Liver",
      "Lymph",
      "Myeloid",
      "Ovary",
      "Panc",
      "Prost",
      "Skin",
      "Stomach"
    ),
    values = circle
  ) +
  labs(x = "tSNE_1", y = "tSNE_2") +
  ggtitle(title)+
  theme(legend.position = 'none')
ggsave(
  file = file,
  plot = g,
  dpi = 500,
  width = 16,
  height = 9
)


# LDAのほう ------------------------------------------------------------------

fileLDA = sprintf(
  "~/Genome/presentation/lda_tsne_%spresen.png",
feature
)
titleLDA = sprintf(
  "LDA featuring : %s",
  feature
)
df2 <- data.frame(matrix(rep(NA, 3), nrow = 1950))[numeric(0), ]
df2 = as.data.frame(cbind(as.factor(type), LDAtsne$Y[, 1] , LDAtsne$Y[, 2]))
df2[, 1] = as.factor(type)
df2[, 4] = as.factor(shape_num)
colnames(df2) <- c("gene_type", "tSNE_1", "tSNE_2", "pch")
g <-
  ggplot(df,
         aes(
           x = df2$tSNE_1,
           y = df2$tSNE_2,
           color = df2$gene_type,
           shape = df2$pch
         )) +
  geom_point() +
  scale_color_manual(
    name = "Cancer detail types",
    labels = sort(unique(type)),
    values = graycolor
  ) +
  scale_shape_manual(
    name = "Cancer types",
    labels = c(
      "Biliary",
      "Bone",
      "Breast",
      "CNS",
      "Eso",
      "Head",
      "Kidney",
      "Liver",
      "Lymph",
      "Myeloid",
      "Ovary",
      "Panc",
      "Prost",
      "Skin",
      "Stomach"
    ),
    values = circle
  ) +
  labs(x = "tSNE_1", y = "tSNE_2") +
  ggtitle(titleLDA) +
  theme(legend.position = 'none')
ggsave(
  file = fileLDA,
  plot = g,
  dpi = 500,
  width = 16,
  height = 9
)
