rm(list=ls())
library(ggplot2)
library(stringr)
filename="gladiator.selected.data.txt"
library(ggplot2)
library("ggpubr")
A<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A)
B=A[1:d[1], 2:d[2]]
#ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
#ZZ=0.000001
ZZ=0
C=B+ZZ

Cname=colnames(A)[2:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
Treat=rep("NA", Clen)
Sid=rep("NA", Clen)
Visit=rep("NA", Clen)
TreatVisit=rep("NA", Clen)
ASFS=rep("NA", Clen)
#age=rep("NA", Clen)
#id2=rep("NA", Clen)
splitname<-strsplit(Cname, "[.]")
for(mm in  1:Clen ){
  Treat[mm]=splitname[[mm]][1]
  if(splitname[[mm]][1] =="A"){
    Treat[mm]="A:ZPT"
  }else if(splitname[[mm]][1] =="B"){
    Treat[mm]="B:PO"
  }else if(splitname[[mm]][1] =="C"){
    Treat[mm]="C:Ctl"
  }
  Visit[mm]=splitname[[mm]][2]
  Sid[mm]=splitname[[mm]][3]
  ASFS[mm]=as.numeric(splitname[[mm]][4])
  #age[mm]=splitname[[mm]][6]
  #id2[mm]=splitname[[mm]][7]
  TreatVisit[mm]=paste0(Treat[mm], Visit[mm])
}
ASFS=as.numeric(ASFS[ASFS !="NA"])

library("plotrix")
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se= std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

for (i in 1:d[1]){
  
  genename=A[i,1]
  gene<-as.numeric(C[i,])
  splitG<-strsplit(as.character(genename), "[.]")
  LLL=length(splitG[[1]])
  genus=splitG[[1]][LLL]
  ccc=as.numeric(min(gene[gene>0&!is.na(gene)]))/10
  RA=gene+ccc
  logRA=log(RA)
  mydata=data.frame(RA,logRA, Treat, Visit, Sid, ASFS, TreatVisit)
  # colnames(mydata)=c("RelativeAbundance", "LogRelativeAbundance", "Treat", "Visit", "SubjectID", "ASFS", "Id2", "TreatVisit")
  df2 <- data_summary(mydata, varname="RA", 
                      groupnames=c("Treat", "Visit"))
  
  
  
  
  png(filename=paste0(i,".",genus, ".png"),  width=600, height=480, res=180)
  p<- ggplot(df2, aes(x=Visit, y=RA,group=Treat, color=Treat)) + 
    geom_line() +
    geom_point(size=3)+
    geom_errorbar(aes(ymin=RA-se, ymax=RA+se), width=.2,
                  position=position_dodge(0.05)) +
    facet_wrap(~Treat)+
    #labs(x="Visit", y = paste0("Abundance: ", genus))+
    labs(x="Visit", y =  paste(genus, "(%)"))+
    theme_classic() +
    theme(legend.position="none")
  
  print(p)
  dev.off()
  
  png(filename=paste0(i,".",genus, ".pair.png"),  width=600, height=480, res=180)
  p1=ggplot(data=mydata,aes(x=Visit,y=RA,color=Treat, group=Sid)) +ylab(paste(genus, "Abundance"))+geom_point(size=1)+geom_line(color="grey")+ theme_classic()+theme(legend.position="none")+facet_wrap(~Treat)
  
  print(p1)
  
  dev.off()
  
}





