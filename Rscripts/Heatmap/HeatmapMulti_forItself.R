#If not already installed
#install.packages("gplots")
#install.packages("devtools")

#Load necessary packages
library("gplots")
library("devtools")
library("RColorBrewer")

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

#Load data
setwd("~/Downloads/nubchi")
fn = "genes.90.addDEG.virulence.heatmap.xls"
data <- read.table(fn, sep = "\t", header = TRUE, row.names = 1)
data <- as.matrix(data)
dim(data)
head(data)
rownames(data)
colnames(data)


#Create fake color side bars
##rownames(data)
#drugclass_colors=sample(c("darkorchid","darkred"), length(drug_names), replace = TRUE, prob = NULL)
#drugcategory_colors=sample(c("green","darkgreen"), length(drug_names), replace = TRUE, prob = NULL)
#rlab=t(cbind(drugclass_colors,drugcategory_colors))
#rownames(rlab)=c("Class","Category")
##colnames(data)
#subtype_colors=sample(c("red","blue","cyan","pink","yellow","green"), length(patient_ids), replace = TRUE, prob = NULL)
#Mcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
#Ncolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
#Tcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
#HER2colors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
#PRcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
#ERcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
#clab=cbind(subtype_colors,Mcolors,Ncolors,Tcolors,HER2colors,PRcolors,ERcolors)
#colnames(clab)=c("Subtype","M","N","T","HER2","PR","ER")
colnames(data)
timeColors=c("lightpink1","lightpink2","lightpink3","lightpink1","lightpink2","lightpink3","lightpink1","lightpink2","lightpink3")
treatColors=c("darkorchid1","darkorchid1","darkorchid1","darkorchid3","darkorchid3","darkorchid3","darkorchid1","darkorchid1","darkorchid1")
virulenceColors=c('gold2','gold2','gold2','gold3','gold3','gold3','gold4','gold4','gold4')
clab=cbind(timeColors,treatColors,virulenceColors)
colnames(clab)=c("TimeSeries Comparison","Treatment Comparison","Virulence Comparison")
colnames(clab)

### make heat-map roughly
#dist2 <- function(x, ...)
#  as.dist(1-cor(t(x), method="pearson"))
myPalette <- colorRampPalette(rev(brewer.pal(5, "RdBu")))
distance.method.list <- c("euclidean", "maximum", "manhattan", "minkowski")
cluster.method.list <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
for(i in 1:length(cluster.method.list)){
  for(j in 1:length(distance.method.list)){
      png(paste(c(cluster.method.list[i],"_", distance.method.list[j],".png"),collapse=""), width=2000, height=2000)
      heatmap.2(log2(data+0.1), col=myPalette, scale="row", trace="none", cexCol=3, cexRow=2, lhei=c(2,20), 
                margins=c(20,20), labRow="",
                hclustfun = function(x) hclust(x,method = cluster.method.list[i]), 
                distfun = function(x) dist(x,method =distance.method.list[j]))
      dev.off()
    }
}

### make heatmap.3 clearly with categories (ColSideColors=clab, RowSideColors=rlab)
main_title="DEG of Virulence pairs"
mydist=function(c) {dist(c,method="minkowski")}
myclust=function(c) {hclust(c,method="complete")}
par(cex.main=1)
heatmap.3(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", dendrogram="both", margins=c(30,10),
          Rowv=TRUE, Colv=TRUE, ColSideColors=clab, symbreaks=FALSE, key=TRUE, symkey=FALSE,
          density.info="none", trace="none", main=main_title, labCol=colnames(data), labRow=FALSE, cexRow=1, col=myPalette,
          ColSideColorsSize=3, RowSideColorsSize=2, KeyValueName=FALSE)
legend("topright",
       legend=c("Highly virulence","Low virulence","","Control","virus","","1 day","3 day","1 week"),
       fill=c("black","grey","white","pink","cyan","white","red","blue","green"), 
       border=TRUE, bty="n", y.intersp = 0.9, cex=0.9)

#legend=c("Basal","LumA","LumB","Her2","Claudin","Normal","","Positive","Negative","NA","","Targeted","Chemo","","Approved","Experimental"),
#fill=c("red","blue","cyan","pink","yellow","green","white","black","white","grey","white","darkorchid","darkred","white","green","darkgreen"), 