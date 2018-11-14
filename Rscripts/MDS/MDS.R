library(Cairo)

raw.exprs.path <- "mature_miRNA_expression.xls"

#1) group 3개로 : WT 1, 2, 3 + CT 1, 9  / WT_ST (TBD180735) / KO_ST (TBD180735)
input_samples <- "WT1:WT2:WT3:CT1:CT9:WT_ST1:WT_ST2:WT_ST3:WT_ST5:KO_ST1:KO_ST2:KO_ST5:KO_ST6"
input_groups <- "1:1:1:1:1:2:2:2:2:3:3:3:3"
MDS_png_file <- "MDS_Group3.png"

#2) group 4개로 : WT 1, 2, 3 + CT 1, 9 / KO (TBD160209) / WT_ST (TBD180735) / KO_ST (TBD180735)
input_samples <- "WT1:WT2:WT3:CT1:CT9:KO1:KO2:KO4:KO5:WT_ST1:WT_ST2:WT_ST3:WT_ST5:KO_ST1:KO_ST2:KO_ST5:KO_ST6"
input_groups <- "1:1:1:1:1:2:2:2:2:3:3:3:3:4:4:4:4"
MDS_png_file <- "MDS_Group4.png"

raw_count = read.table(file=raw.exprs.path, header=TRUE, row.names=1, sep='\t', quote = "")
#colnames(raw_count)
samples <- make.names(unlist(strsplit(input_samples,":")))
sel_col <- samples
groups <- unlist(strsplit(input_groups,":"))
group <- groups

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


# MDS plot with expression data
t_data <- t(raw_count[,sel_col])
d <-dist(log(t_data+0.00001,2),method='euclidean')
fit <- cmdscale(d,eig=TRUE,k=2)
x <- fit$points[,1]
y <- fit$points[,2]
CairoPNG(MDS_png_file, width=800, height=800)
plot(x,y,xlab="Coordinate 1", ylab="Coordinate 2", main="MDS(log2(fpkm))", type="n")
text(x,y,labels = row.names(t_data), cex=.9, col=string.to.color(group))
dev.off()
