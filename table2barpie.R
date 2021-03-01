############################################
#need to try split panel
#########################################
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outname <- args[2]
#filename="5945_Short_Genus_Mean"
#outname="Mcup_Genus"
A<-read.table(filename,header=TRUE, sep="\t",stringsAsFactors=F)
d<-dim(A)
newsum=A[,2:d[2]]
bugs=A[,1]
row.names(newsum)=A[,1]
ttt=colnames(A)[1]
ddd=dim(newsum)
library(RColorBrewer)
col1=brewer.pal(9, "Set1")
col2=brewer.pal(12, "Set3")
col3=brewer.pal(9, "Pastel1")
col4=brewer.pal(8, "Accent")
col5=brewer.pal(8, "Dark2")

cols=c(col1, col2, col3, col4, col5)[1:ddd[1]]
#cols=colors()[1:ddd[1]]
#cols=colors(c(2,3,8,15,18,19,24,26,31,32,33,37,42,43,47,52,62,65,68,73,74,76,81,82,83,84, 86,96,109,116,121,136,142,384,387,388,400))

for(i in 1:ddd[2]){
      B=newsum[,i]
      percentlabels<- round(B, 2)
#      pielabels<- paste(rownames(newsum),":", percentlabels, "%", sep="")
      pielabels<- paste(percentlabels, "%", sep="")
      pielabels[percentlabels<1]=""
      png(filename=paste0(outname, ".",colnames(newsum)[i],".pie.png"), width=1600, height=800,res=120 )      
      pie(B, main=colnames(newsum)[i],col=cols,  labels=pielabels, cex=0.7)
       dev.off()
       png(filename=paste0(outname, ".",colnames(newsum)[i],".pie.nolabel.png"), width=600, height=600,res=120 )      
      pie(B, main=colnames(newsum)[i],col=cols, labels=NA)
       dev.off()

}


#png(filename=paste0(outname,".bar.png"), width=300, height=800,res=120 ) 

png(filename=paste0(outname,".bar.png"), width=1200, height=800,res=120 ) 
par(mar=c(9.8,4.1,4.1,2.1))
barplot(as.matrix(newsum), col=cols, width=2, ylab="relative abundance", las=2, cex.names=1
)
dev.off()

png(filename=paste0(outname,".lengend.png"), width=800, height=1200,res=120)
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
legend("topleft",bugs, cex=0.8, fill=cols, title=ttt)
dev.off()

