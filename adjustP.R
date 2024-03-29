args <- commandArgs(trailingOnly = TRUE)
#print(args)
filename <- args[1]

#filename="gene.community.relab10.stat.xls"
A<-read.table(filename,header=TRUE, sep="\t",stringsAsFactors=F)
d=dim(A)
cols=c(3,4,5,6,7,8,9,10,11,12,13)
Ulen=length(cols)
df<-data.frame(matrix(ncol=Ulen, nrow=d[1]))
colnames(df)=paste0("fdr.", colnames(A)[cols])
rownames(df)=A[,1]
for(i in 1:Ulen){
 colN=cols[i]
 p=A[,colN]
 df[,i]=p.adjust(p, method="fdr")
}
write.table(df, file=paste0(filename, ".adjustp"), sep = "\t",eol = "\n", na = "NA", row.names = TRUE,col.names = TRUE)
