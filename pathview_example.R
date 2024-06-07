# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview")


library(pathview)

mypathway <-"00020"
genes <-c("K01596", "K00169", "K00627", "K01681")
logFC<-rnorm(length(genes),-0.5,1)
names(logFC)<-genes
pathview(gene.data=logFC,species="ko",pathway=mypathway)


mypathway<-"vvi03060"
###This example use entrenz ID" 
 genes<-c("100241050","100243802","100244217","100244265","100247624","100247887","100248517","
 100248990","100250268","100250385","100250458","100251379","100252350","100252527","100252725","
 100252902","100253826","100254350","100254429","100254996","100255515","100256046","100256113","
 100256412","100256941","100257568","100257730","100258179","100258854","100259285","100259443","
 100260422","100260431","100261219","100262919","100263033","100264739","100265371","100266802","
 100267343","100267692","100852861","100853033","100854110","100854416","100854647","100855182","
 104879783","109122671")
logFC<-rnorm(length(genes),-0.5,1)
names(logFC)<-genes
pathview(gene.data=logFC,species="vvi",pathway=mypathway)
