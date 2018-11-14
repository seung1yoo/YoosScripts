
library(Cairo)

string.to.color = function(string, colors=NULL){
  if (!is.null(colors)){
    if (length(colors)!=length(unique(string))){
      break("The number of colors must be equal to the number of unique elements.")
    }
    else {
      conv = cbind(unique(string), colors)
    }
  } else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN=function(x){conv[which(conv[,1]==x),2]}))
}

stdin <- commandArgs(TRUE)

raw.exprs.path <- stdin[1]
output.dir <- stdin[2]
label <- stdin[3]
samples <- make.names(unlist(strsplit(stdin[4],":")))
groups <- unlist(strsplit(stdin[5],":"))
if(length(stdin) > 5){
  pairs <- unlist(strsplit(stdin[6],":"))
}

#raw.exprs.path <- '/BiO/BioProjects/ISH-Human-SmallRNAseq-2015-12/Result/report/expression/mature_miRNA_expression.xls'
#output.dir <- '/BiO/BioProjects/ISH-Human-SmallRNAseq-2015-12/Result/report/differential_expression'
#label <- 'DEG_001'
#samples <- make.names(unlist(strsplit("HCC 02-03-2015 PVTT:TISV NT 20150924 2-1",":")))
#groups <- unlist(strsplit("T:N",":"))

type <- 'edger'
#type <- 'deseq'

type
samples

MA_png_file = paste(output.dir,"/",label,".MA.png",sep="")
norm_count_xls = paste(output.dir,"/",label,".norm_count.xls",sep="")
deg_list_xls = paste(output.dir,"/",label,".DEG.xls",sep="")
MDS_png_file = paste(output.dir,"/",label,".MDS.png",sep="")

sel_col <- samples
group <- groups

raw_count=read.table(file=raw.exprs.path,header=TRUE, row.names=1, sep='\t', quote = "")

#Constructing TCC class object
library(TCC)

if (exists("pair")){
  cl <- data.frame(group = group, pair = pair)
}else{
  cl <- data.frame(group = group)
}
tcc <- new ("TCC", raw_count[,sel_col], cl) # use specific samples
#tcc <- new ("TCC",raw_count,group) # use all samples

#Filtering low-count genes
dim(tcc$count)
tcc <- filterLowCountGenes(tcc)
dim(tcc$count)
library(DESeq)

# DEGES/edge-R
if (exists("pair") & type == "edger"){
  tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger",
                         iteration=1, FDR=0.1, floorPDEG=0.05, paired = TRUE)
  tcc <- estimateDE(tcc, test.method="edger", FDR=0.1, paired = TRUE)
}else if (type == "edger"){
  tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger",
                         #iteration=1, FDR=0.1, floorPDEG=0.05)
  #tcc <- calcNormFactors(tcc, norm.method="deseq2", test.method="deseq2",
                         iteration=1, FDR=0.1, floorPDEG=0.05)
  tcc <- estimateDE(tcc, test.method="edger", FDR=0.1)
  #tcc <- estimateDE(tcc, test.method="deseq2", FDR=0.1)
}

# DEGES/DESeq
if (exists("pair") & type == "deseq"){
  tcc <- calcNormFactors(tcc, norm.method="deseq", test.method="deseq",
                         iteration=1, FDR=0.1, floorPDEG=0.05, paired = TRUE)
  tcc <- estimateDE(tcc, test.method="deseq", FDR=0.1, paired = TRUE)
}else if (type == "deseq"){
  tcc <- calcNormFactors(tcc, norm.method="deseq", test.method="deseq",
                         iteration=1, FDR=0.1, floorPDEG=0.05)
  tcc <- estimateDE(tcc, test.method="deseq", FDR=0.1)
}

table(tcc$estimatedDEG)

# The plot function generates an M-A plot
CairoPNG(MA_png_file, width=640, height=480)
plot(tcc, median.lines = TRUE, cex=0.4)
dev.off()

# The getNormalizedData function can be applied to the TCC class object after the normalization factors have been calculated
normalized.count <- getNormalizedData(tcc)
write.table(round(normalized.count,3), norm_count_xls,
            sep="\t",quote=F,col.names=NA, row.names=T)

# The summary statistics for top-ranked genes
result <- getResult(tcc, sort=TRUE)
write.table(result, deg_list_xls,
            sep="\t",quote=F,col.names=T,row.names=F)

#library(ggplot2)

#An MA-plot is a plot of log-fold change (M-values, i.e. the log of the ratio of level counts for each gene between two samples) against the log-average (A-values, i.e. the average level counts for each gene across the two samples).
#The MA-plot is a useful to visualize reproducibility between samples of an experiment. From a MA-plot one can see if normalization is needed.
#In MA plot, genes with similar expression levels in two samples will appear around the horizontal line y = 0. A lowess fit (in red) is plotted underlying a possible trend in the bias related to the mean expression.
#A = result$a.value
#M = result$m.value
#df = data.frame(A,M)
#ggplot(df, aes(x = A, y = M)) + geom_point(size = 1.5, alpha = 1/5) +
#  geom_hline(color = "blue3") + stat_smooth(se = FALSE, method = "loess", color = "red3")

# MDS plot with expression data
t_data <- t(raw_count[,sel_col])
d <-dist(log(t_data+0.00001,2),method='euclidean')
fit <- cmdscale(d,eig=TRUE,k=2)
x <- fit$points[,1]
y <- fit$points[,2]
CairoPNG(MDS_png_file, width=800, height=800)
plot(x,y,xlab="Coordinate 1", ylab="Coordinate 2", main="MDS(log2(fpkm))", type="n")
text(x,y,labels = row.names(t_data), cex=.7, col=string.to.color(group))
dev.off()
