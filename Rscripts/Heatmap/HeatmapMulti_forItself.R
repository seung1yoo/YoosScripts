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
fn = "genes.90.addDEG.size.Int.Only.Heatmap.EXP.AVG.xls"
data <- read.table(fn, sep = "@", header = TRUE, row.names = 1)
data <- read.csv(fn, sep = ",", row.names = 1)
data <- as.matrix(data)
dim(data)
head(data)
rownames(data)
colnames(data)

#column names change
colnames(data)<-c("20C_virus_1day","20C_virus_3day","20C_virus_7day","20C_control_1day","20C_control_3day","20C_control_7day","13C_virus_1day",
                  "13C_virus_3day","13C_virus_7day","13C_control_1day","13C_control_3day","13C_control_7day")

#extracting specific columns from a data frame
colnames(data)
colnames(data[,c(1,2,3,7,8,9)])
data <- data[,c(1,2,3,7,8,9)]
colnames(data)
head(data,1)

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
timeColors=c("lightpink1","lightpink2","lightpink3","lightpink1","lightpink2","lightpink3","lightpink1","lightpink2","lightpink3","lightpink1","lightpink2","lightpink3")
treatColors=c("darkorchid","darkorchid","darkorchid","darkorchid1","darkorchid1","darkorchid1","darkorchid","darkorchid","darkorchid","darkorchid1","darkorchid1","darkorchid1")
virulenceColors=c('gold2','gold2','gold2','gold2','gold2','gold2','gold4','gold4','gold4','gold4','gold4','gold4')
clab=cbind(timeColors,treatColors,virulenceColors)
clab
colnames(clab)=c("TimeSeries","Treatment","Size")
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
                margins=c(40,20), labRow="",
                hclustfun = function(x) hclust(x,method = cluster.method.list[i]), 
                distfun = function(x) dist(x,method =distance.method.list[j]))
      dev.off()
    }
}

### make heatmap.3 clearly with categories (ColSideColors=clab, RowSideColors=rlab)
main_title=""
mydist=function(c) {dist(c,method="minkowski")}
myclust=function(c) {hclust(c,method="mcquitty")}
myPalette <- colorRampPalette(rev(brewer.pal(5, "RdBu")))
par(cex.main=1)
heatmap.3(log2(data+0.1), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", dendrogram="both", margins=c(15,10),
          Rowv=TRUE, Colv=TRUE, ColSideColors=clab, symbreaks=FALSE, 
          key=TRUE, keysize = 1, symkey=FALSE,
          density.info="none", trace="none", main=main_title, labCol=colnames(data), labRow=FALSE, cexRow=1, cexCol=1.3, 
          col=myPalette, ColSideColorsSize=3, RowSideColorsSize=2, KeyValueName=FALSE)
# legend : Specify legend position by keywords
# http://www.sthda.com/english/wiki/add-legends-to-plots-in-r-software-the-easiest-way
legend("topright", inset = .02,
       legend=c("20C","13C","","Control","virus","","1d","3d","7d"),
       fill=c("gold2","gold4","white","darkorchid1","darkorchid","white","lightpink1","lightpink2","lightpink3"), 
       border=FALSE, bty="n", angle = 90, y.intersp = 1, cex=0.7)

#legend=c("Basal","LumA","LumB","Her2","Claudin","Normal","","Positive","Negative","NA","","Targeted","Chemo","","Approved","Experimental"),
#fill=c("red","blue","cyan","pink","yellow","green","white","black","white","grey","white","darkorchid","darkred","white","green","darkgreen"), 