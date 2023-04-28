rm(list=ls())
library(ggplot2)
library(stringr)
filename="GO_Select_data.xls"
library(ggplot2)
library("ggpubr")
A<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A)
B=A[1:d[1], 4:d[2]]
#ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
#ZZ=0.000001
ZZ=0
C=B+ZZ

Cname=colnames(A)[4:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
Treat=rep("NA", Clen)
#Sid=rep("NA", Clen)
#Visit=rep("NA", Clen)
#TreatVisit=rep("NA", Clen)
#ASFS=rep("NA", Clen)
#age=rep("NA", Clen)
#id2=rep("NA", Clen)
GroupName=rep("NA", Clen)
splitname<-strsplit(Cname, "[.]")
for(mm in  1:Clen ){
  Treat[mm]=splitname[[mm]][2]
  if((splitname[[mm]][2] =="1ND") ||
    (splitname[[mm]][2] =="2MID") || 
    (splitname[[mm]][2] =="3D")){
    GroupName[mm]="Dandruff_Etiology"
  }else if((splitname[[mm]][2] =="BBL")||(splitname[[mm]][2] =="BW3")){
    GroupName[mm]="PO_Shampoo"
  }else if((splitname[[mm]][2] =="CBL")||(splitname[[mm]][2] =="CW3")){
    GroupName[mm]="Control_Shampoo"
  }
  #Visit[mm]=splitname[[mm]][2]
  #Sid[mm]=splitname[[mm]][3]
  #ASFS[mm]=as.numeric(splitname[[mm]][4])
  #age[mm]=splitname[[mm]][6]
  #id2[mm]=splitname[[mm]][7]
  #TreatVisit[mm]=paste0(Treat[mm], Visit[mm])
}
#ASFS=as.numeric(ASFS[ASFS !="NA"])
Treat[Treat=="1ND"]="1NoDandruff"
Treat[Treat=="3D"]="3Dandruff"
Treat[Treat=="2MID"]="2Intermediate"
Treat[Treat=="BBL"]="4PO_Baseline"
Treat[Treat=="BW3"]="5PO_Week3"
Treat[Treat=="CBL"]="6Control_Baseline"
Treat[Treat=="CW3"]="7Control_Week3" 

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
  goid=gsub(":", "_", A[i,1])
  goid=gsub(" ", "_", goid)
  gene<-as.numeric(C[i,])
  #splitG<-strsplit(as.character(genename), "[.]")
  #LLL=length(splitG[[1]])
  #genus=splitG[[1]][LLL]
  ccc=as.numeric(min(gene[gene>0&!is.na(gene)]))/10
  RA=gene+ccc
  logRA=log(RA)
 
  mydata=data.frame(RA,logRA, Treat, GroupName)
  # colnames(mydata)=c("RelativeAbundance", "LogRelativeAbundance", "Treat", "Visit", "SubjectID", "ASFS", "Id2", "TreatVisit")
  df2 <- data_summary(mydata, varname="RA", 
                      groupnames=c("Treat", "GroupName"))
  
 
  
  
  png(filename=paste0(i,".",substr(goid, 1,80), ".png"),  width=1380, height=480, res=180)
  p<- ggplot(df2, aes(x=Treat, y=RA,group=GroupName, color=GroupName)) + 
    geom_line() +
    geom_point(size=5)+
    geom_errorbar(aes(ymin=RA-se, ymax=RA+se), width=.2,
                  position=position_dodge(0.05)) +
    labs(x="", y ="Relative Abundance")+
    theme_classic() +
    theme(legend.position="none")+ggtitle(genename)

  
  
  print(p)
  dev.off()
  
  #png(filename=paste0(i,".",genus, ".pair.png"),  width=600, height=480, res=180)
  #p1=ggplot(data=mydata,aes(x=Visit,y=RA,color=Treat, group=Sid)) +ylab(paste(genus, "Abundance"))+geom_point(size=1)+geom_line(color="grey")+ theme_classic()+theme(legend.position="none")+facet_wrap(~Treat)
  
  #print(p1)
  
  #dev.off()
  
}





