d = fread(
  "mutation_small.csv",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  sep = ","
)


data = subset(data, data$FILTER %in% "PASS")
library("rospca")

X <- dataGen(m=1, n=100, p=10, eps=0.2, bLength=4)$data[[1]]

resR <- robpca(X, k=2)
biplot(resR,)
print(
  ggbiplot(
    resR,
    obs.scale = 1,
    var.scale = 1,
    # groups = label_small,
    ellipse = TRUE,
    circle = TRUE,
    # color = label_small
  )
  + ggtitle("title")
)

diagPlot(resR)



install.packages("cvreg_0.2.0.tar.gz",repos = NULL, type = "source")


