##########################
rm(list=ls())
library(stringr)
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

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
#filename="metaphlan.relab10.7"
#filename="MulSanne.species10.txt"
library(ggplot2)
library("ggpubr")
my.t.test.p.value <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}

my.wilcox.p.value <- function(...) {
    obj<-try(wilcox.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}
truefc<-function(VVV){
	XXX=VVV
	if(VVV==0){
	    XXX=NA
   	}else if(VVV<1){
	    XXX=-1/VVV
    	}
	return(XXX)
}

A<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A);
B=A[1:d[1], 2:d[2]]
ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
C=B+ZZ

Cname=colnames(A)[2:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "_")


batch=rep("NA", Clen)
Snumber=rep("NA", Clen)
SID=rep("NA", Clen)
Visit=rep("NA", Clen)
Gender=rep("NA", Clen)
Disease=rep("NA", Clen)
for(mm in  1:Clen ){
  Disease[mm]=splitname[[mm]][2]
  Visit[mm]=splitname[[mm]][3]
  SID[mm]=splitname[[mm]][4]
  Gender[mm]=splitname[[mm]][5]
  batch[mm]=splitname[[mm]][6]
  Snumber[mm]=splitname[[mm]][7]
}

DiseaseVisit=paste0(Disease, "_", Visit)
DiseaseVisitSID=paste0(Disease, "_", Visit,  "_", SID)
custom_order <- c(
  "HEALTHY_BL", "HEALTHY_W4","HEALTHY_W8", "UNHEALTHY_BL", "UNHEALTHY_W4",    "UNHEALTHY_W8"    
  
)
DiseaseVisit <- factor(DiseaseVisit, levels = custom_order)

print(paste("genename","genus2","percentZeroNA","KP_Disease"," KP_SID","KP_Visit","KP_Gender","KP_DiseaseVisit","Pwilcox_UBL_HBL","Pwilcox_UW4_UBL","Pwilcox_UW8_UBL"," truefc_UBL_HBL","truefc_UW4_UBL","truefc_UW8_UBL","MUBL","MHBL","MUW4","MHW4","MUW8","MHW8",   sep=","))
#
for (i in 1:d[1]){
    genename=A[i, 1]
    if(is.na(genename)){next;}
    gene<-as.numeric(C[i,])
    splitG<-strsplit(as.character(genename), "[.|]")
    LLL=length(splitG[[1]])
    genus=splitG[[1]][LLL]
    genus1=str_replace(genus, "[kposgfc]__", "")
    genus2=str_replace(genus1, "_", " ")
    mydata=data.frame(gene,log(gene), Visit,Disease,SID, DiseaseVisit, DiseaseVisitSID, Gender,batch, Snumber)
    
    KP_Visit=kruskal.test(gene ~ Visit, data = mydata)$p.value
    KP_Disease=kruskal.test(gene ~ Disease, data = mydata)$p.value
    KP_DiseaseVisit=kruskal.test(gene ~ DiseaseVisit, data = mydata)$p.value
    KP_SID=kruskal.test(gene ~ SID, data = mydata)$p.value
    KP_Gender=kruskal.test(gene ~ Gender, data = mydata)$p.value
    H_BL<-gene[DiseaseVisit=="HEALTHY_BL"]; MHBL=mean(H_BL)
    H_W4<-gene[DiseaseVisit=="HEALTHY_W4"]; MHW4=mean(H_W4)
    H_W8<-gene[DiseaseVisit=="HEALTHY_W8"]; MHW8=mean(H_W8)
    U_BL<-gene[DiseaseVisit=="UNHEALTHY_BL"]; MUBL=mean(U_BL)
    U_W4<-gene[DiseaseVisit=="UNHEALTHY_W4"]; MUW4=mean(U_W4)
    U_W8<-gene[DiseaseVisit=="UNHEALTHY_W8"]; MUW8=mean(U_W8)
    
  
    truefc_UBL_HBL=truefc(MUBL/MHBL)
    truefc_UW4_UBL=truefc(MUW4/MUBL)
    truefc_UW8_UBL=truefc(MUW8/MUBL)
    Pwilcox_UBL_HBL<-my.wilcox.p.value(U_BL, H_BL, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_UW4_UBL<-my.wilcox.p.value(U_W4, U_BL, na.rm=TRUE,paired = TRUE, alternative = "two.sided")
    Pwilcox_UW8_UBL<-my.wilcox.p.value(U_W8, U_BL, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
  
    percentZeroNA=(sum(B[i,]==0)+sum(is.na(B[i,])))/(d[2]-1)
    
    print(paste(genename,genus2,percentZeroNA, KP_Disease,  KP_SID, KP_Visit, KP_Gender, KP_DiseaseVisit, 
                Pwilcox_UBL_HBL, Pwilcox_UW4_UBL, Pwilcox_UW8_UBL,  truefc_UBL_HBL, truefc_UW4_UBL, truefc_UW8_UBL, 
                MUBL, MHBL, MUW4, MHW4, MUW8, MHW8,   sep=","))
    pvalues <- c(Pwilcox_UBL_HBL, Pwilcox_UW4_UBL, Pwilcox_UW8_UBL, KP_DiseaseVisit)
    formatted_pvalues <- sprintf("%.2f", pvalues)
     if(percentZeroNA<0.75){
      if(min(formatted_pvalues)<=0.05){
        
        df2 <- data_summary(mydata, varname = "gene", groupnames = c("DiseaseVisit", "Disease", "Visit"))
        
        # Use factor to set custom order
        df2$DiseaseVisit <- factor(df2$DiseaseVisit, levels = custom_order)
        
        p <- ggplot(df2, aes(y = gene, x = DiseaseVisit, group = Disease, color = Disease)) +
          geom_point(size = 5) + geom_line()+
          geom_errorbar(aes(ymin = gene - se, ymax = gene + se), width = .2,
                        position = position_dodge(0.05)) +
          labs(y = paste(genename, "relative abundance"), x = "") +
          theme_classic() +
          theme(legend.position = "none") +
          ggtitle(paste(
            "pUBL_HBL:", formatted_pvalues[1],
            ";pUW4_UBL:", formatted_pvalues[2],
            ";pUW8_UBL:", formatted_pvalues[3],
            ";KpDiseaseVisit:", formatted_pvalues[4],
            "; ", genename
          )) +
          theme(plot.margin = margin(0.5, 0.5, 0.5, 1.5, "cm"), plot.title = element_text(size = 10))  # Adjust the title font size
        
        # Save the plot
        ggsave(paste0(genename, ".png"), p, width = 7.5, height = 4)
        
       
      }
    }
    
}
