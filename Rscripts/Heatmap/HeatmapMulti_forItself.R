setwd("/Volumes/TBI_siyoo/TBI_Research/04_GNU/01.KimEuiYeon_TBD160336/5.GSEAv2/goInfo/CHLOROPHYLLIDE_A_OXYGENASE_[OVERALL]_ACTIVITY")
fn = "exp.txt"
setwd("/Volumes/TBI_siyoo/TBI_Research/04_GNU/02.BackDongWon_TBD160282/2.Analysis/report/GSEAv2/Ara_gsea_GO/Top5GO/")
fn = "total_exp.redupl.txt"

setwd("/Volumes/TBI_siyoo/TBI_Research/DKU_HGD_Tapes_TBD150282/RNAseq/DKU_Han-Venerupis-2016-03_V1_Ref_P005/GSEA/Tapes_gsea_v2/heatmap/")
fn = "Tapes.exp.txt"
data <- data[c("AGAMI","KWANJA","BAL")]

data <- read.table(fn, sep="\t", header = TRUE, row.names = 1, quote = "")
data <- as.matrix(data)
dim(data)
head(data)
library(gplots)
library("RColorBrewer")

dist2 <- function(x, ...)
  as.dist(1-cor(t(x), method="pearson"))

myPalette <- colorRampPalette(rev(brewer.pal(5, "RdBu")))
cluster.method.list <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
distance.method.list <- c("euclidean", "maximum", "manhattan", "minkowski")

for(i in 1:length(cluster.method.list)){
  for(j in 1:length(distance.method.list)){
      png(paste(c(cluster.method.list[i],"_", distance.method.list[j],".png"),collapse=""), width=2000, height=2000)
      #heatmap.2(log2(data+0.1), col=myPalette, labRow = "", scale="row", trace="none", key=FALSE, cexCol=3, cexRow=2, lhei=c(2,20), margins=c(20,20),   hclustfun = function(x) hclust(x,method = cluster.method.list[i]), distfun = function(x) dist(x,method =distance.method.list[j]))
      heatmap.2(log2(data+0.1), col=myPalette, scale="row", trace="none", cexCol=3, cexRow=2, lhei=c(2,20), 
                margins=c(20,20), labRow="",
                hclustfun = function(x) hclust(x,method = cluster.method.list[i]), 
                distfun = function(x) dist(x,method =distance.method.list[j]))
      dev.off()
    }
}

png("heatmap.png",width=1024,height=1024)
heatmap.2(log2(data+0.1), col=myPalette, trace="none", cexRow=0.7, cexCol=1, margins=c(7,8), labRow="")
dev.off()
