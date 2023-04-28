rm(list=ls())
library(ggplot2)

filename="vox1.top7.data.txt"
outname="Vox1.metaphlan4.7"
A<-read.table(filename,header=TRUE, sep="\t",stringsAsFactors=F)
d<-dim(A)
newsum=A[,8:d[2]]
avesum=A[,2:6]
allsum=A[,2:d[2]]
bugs=A[,1]
ddd=dim(newsum)
library(RColorBrewer)
col1=brewer.pal(9, "Set1")
col2=brewer.pal(12, "Set3")
col3=brewer.pal(9, "Pastel1")
col4=brewer.pal(8, "Accent")
col5=brewer.pal(8, "Dark2")

cols=c(col1, col2, col3, col4, col5)[1:ddd[1]]

png(filename=paste0(outname,".bar.png"), width=1800, height=800,res=120 ) 
par(mar=c(9.8,4.1,4.1,2.1))
barplot(as.matrix(newsum), col=cols, width=2, ylab="relative abundance", las=2, cex.names=0.8, border="NA")+ theme_classic()+ theme(axis.text.x = element_text(size = 1))
dev.off()

png(filename=paste0(outname,".avebar.png"), width=300, height=800,res=120 ) 
par(mar=c(9.8,4.1,4.1,2.1))
barplot(as.matrix(avesum), col=cols, width=2, ylab="relative abundance", las=2, cex.names=0.8, border="NA")+ theme_classic()+ theme(axis.text.x = element_text(size = 1))
dev.off()

png(filename=paste0(outname,".allbar.png"), width=2000, height=800,res=120 ) 
par(mar=c(9.8,4.1,4.1,2.1))
barplot(as.matrix(allsum), col=cols, width=2, ylab="relative abundance %", las=2, cex.names=0.8, border="NA")+ theme_classic()+ theme(axis.text.x = element_text(size = 0.2))
dev.off()

png(filename=paste0(outname,".lengend.png"), width=800, height=1200, res=120)
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
legend("topleft",bugs, cex=0.8, fill=cols, title=colnames(A)[1])
dev.off()


