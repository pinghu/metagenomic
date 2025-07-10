#options(echo=TRUE) # if you want see commands in output file
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
print(args)


filename=args[1]
outname=args[2]
#filename ="Mucositis_merged_data_small.txt"
#outname ="Mucositis-biomarker"
#filename="Mucositis_merged_data_select3.txt"
#outname="genus_biomarker"
#filename ="Mucositis_select2_data.txt"
#outname ="Mucositis"

rm(args)

A<-read.table(filename, sep="\t", header=TRUE)

d <- dim(A);
B <- data.matrix(A[1:d[1],2:d[2]]);
ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
test_data=log10(B+ZZ)
xx<-t(test_data)

#B<-cor(xx)
colnames(xx) =A[,1]
rownames(xx)=colnames(A)[2:d[2]]

library("Hmisc")
rCOR=rcorr(as.matrix(xx), type="spearman")

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
#corrplot(res, type = "upper", order = "hclust", 
#         tl.col = "black", tl.srt = 45)
#tiff(filename = paste0("Fig4a.tiff"), width = 8000, height = 8000, res=300)
#corrplot.mixed(rCOR$r, p.mat=rCOR$P, sig.level=0.05, diag='u', insig = "blank",  order="hclust", upper = "ellipse", lower.col = "black", tl.pos="lt")
#dev.off()

CorResult=cor(as.matrix(xx), method="spearman")
rCOR2=rcorr(as.matrix(xx), type="spearman")
testRes = cor.mtest(as.matrix(xx), method="spearman" , conf.level = 0.95)

tiff(filename = paste0(outname, ".spearmanP.cluster.tiff"), width = 4298, height = 4298, res=300)
corrplot.mixed(corr = CorResult, p.mat = testRes$p, sig.level=0.05,  tl.pos="lt", tl.col = "blue4", lower.col = "black",insig = "blank",font = 3, upper = "square", lower  = "number"
               #, order="hclust"
               )
dev.off()


tiff(filename = paste0(outname, ".spearmanP.cluster.1.tiff"), width = 4298, height = 4298, res=300)
corrplot.mixed(corr = CorResult, p.mat = testRes$p, sig.level=0.1,  tl.pos="lt", tl.col = "blue4", lower.col = "black",insig = "blank",font = 3, upper = "square", lower  = "number"
               #, order="hclust"
)
dev.off()




tiff(filename = paste0(outname, ".spearman.tiff"), width = 4888, height = 4888, res=300)
#corrplot(corr = CorResult, method = "ellipse", type = "upper",
#         col = colorRampPalette(c("#a87963", "#4b61ba"))(10))
#corrplot.mixed(corr = CorResult, lower = "square", upper = "number",
#               lower.col = colorRampPalette(c("#a87963", "#4b61ba"))(10),
#               upper.col = colorRampPalette(c("#a87963", "#4b61ba"))(10))
corrplot.mixed(corr =CorResult,  lower = "square", upper = "number", tl.pos="lt")
dev.off()
#library(metan)
#ALL <-corr_coef(xx)