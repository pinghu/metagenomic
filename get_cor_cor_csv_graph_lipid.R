#options(echo=TRUE) # if you want see commands in output file
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
print(args)
#filename ="RA6.select.SG.clean"
#outname ="RA6SG"
filename=args[1]
outname=args[2]

rm(args)

A<-read.table(filename, sep="\t", header=TRUE)

d <- dim(A);
test_data <- data.matrix(A[1:d[1],2:d[2]]);
xx<-t(test_data)

#B<-cor(xx)
colnames(xx) =A[,1]
rownames(xx)=colnames(A)[2:d[2]]

library("Hmisc")
rCOR=rcorr(as.matrix(xx))
write.table(rCOR$P,file=paste0(outname,'.corP.csv'),sep=",",row.names=A[,1], col.names = A[,1])
write.table(rCOR$r,file=paste0(outname,'.corR.csv'),sep=",",row.names=A[,1], col.names = A[,1])
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

ff<-flattenCorrMatrix(rCOR$r, rCOR$P)
write.table(ff,file=paste0(outname,'.corFlat.tsv'),sep="\t")

library(corrplot)
res <- cor(xx)
round(res, 2)

jpeg(filename = paste0(outname,".Cor.jpeg"), width = 4860, height = 4860, res=300)
corrplot.mixed(rCOR$r, p.mat=rCOR$P, rect.hc = NA,sig.level=0.05, insig = "blank",  order="original",  upper = "ellipse", lower.col = "black", tl.pos="lt")
dev.off()
