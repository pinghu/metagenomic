
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
print(args)
filename <- args[1]
#metafile <-args[2]
outname <- filename
rm(args)
metafile="../data/freshpad.otherdata"
#filename="../data/kraken_fungi18DB_summary.clean.7.recal.Percent.RShort.xls"
scalefile="../data/freshpad.Absolute.ScalFactor"
D<-read.table(scalefile, sep="\t", header=TRUE)
scaleF=as.numeric(D[1,2:dim(D)[2]])

A<-read.table(filename, sep="\t", header=TRUE)
X <- A[,!grepl("Control", colnames(A))]
d<-dim(X)
X2<-X
for (i in 1:d[1]){
  X2[i,2:d[2]]<-as.numeric(X[i,2:d[2]])*scaleF
}
write.table(X2,paste0(filename, ".scaled"),sep="\t",row.names=FALSE)

C<-read.table(metafile, sep="\t", header=TRUE)
colnames(X2) =colnames(C)
B0<-rbind(C, X2)
d<-dim(B0)
B<- data.matrix(B0[1:d[1],2:d[2]])
ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
test_data=log10(B+ZZ)
xx<-t(test_data)
colnames(xx) =B0[,1]
rownames(xx)=colnames(B0)[2:d[2]]

#B<-cor(xx)

library("Hmisc")
rCOR=rcorr(as.matrix(xx), type="spearman")
M=cor(as.matrix(xx), method="spearman")

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
write.table(ff,file=paste0(outname,'.corFlatS.tsv'),sep="\t")


rCORp=rcorr(as.matrix(xx), type="pearson")
Mp=cor(as.matrix(xx), method="pearson")
ff<-flattenCorrMatrix(rCORp$r, rCORp$P)
write.table(ff,file=paste0(outname,'.corFlatP.tsv'),sep="\t")
