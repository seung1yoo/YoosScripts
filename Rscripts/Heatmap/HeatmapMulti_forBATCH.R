
## R CMD BATCH --no-save --no-restore '--arg file outDir' HeatmapMulti_forBATCH.R HeatmapMulti_forBATCH.Rout

args <- commandArgs(TRUE)
fn <- args[1]
print (fn)
outDir <- args[2]
print (fn)

data <- read.table(fn, sep="\t", header = TRUE, row.names = 1, quote = "")
data <- as.matrix(data)
dim(data)
head(data)
library(gplots)
library("RColorBrewer")

dist2 <- function(x, ...)
  as.dist(1-cor(t(x), method="pearson"))

myPalette <- colorRampPalette(rev(brewer.pal(5, "RdBu")))
#cluster.method.list <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
#distance.method.list <- c("euclidean", "maximum", "manhattan", "minkowski")
cluster.method.list <- c("complete")
distance.method.list <- c("euclidean")

for(i in 1:length(cluster.method.list)){
  for(j in 1:length(distance.method.list)){
      png(paste(c(outDir, paste(c(cluster.method.list[i],"_", distance.method.list[j],".png"),collapse="")), collapse="/"), width=2000, height=2000)
      heatmap.2(log2(data+0.1), col=myPalette, trace="none", cexCol=3, cexRow=2, lhei=c(2,20), margins=c(20,20),   hclustfun = function(x) hclust(x,method = cluster.method.list[i]), distfun = function(x) dist(x,method =distance.method.list[j]))
      #heatmap.2(log2(data+0.1), col=myPalette, scale="row", trace="none", cexCol=3, cexRow=2, lhei=c(2,20), margins=c(20,20),   hclustfun = function(x) hclust(x,method = cluster.method.list[i]), distfun = function(x) dist(x,method =distance.method.list[j]))
      dev.off()
    }
}

