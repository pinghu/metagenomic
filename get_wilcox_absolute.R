######################################################
#I also need to add the correlation into the picture
###need to remove the # from the title page
#######################################################
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(stringr)
args <- commandArgs(trailingOnly = TRUE)
#print(args)
filename <- args[1]
#filename="../data/metaphlan.relab10.7"
metafile="../data/freshpad.Absolute.ScalFactor"
D<-read.table(metafile, sep="\t", header=TRUE)
scaleF=as.numeric(D[1,2:dim(D)[2]])
#rm(args)
#filename="Vox1_taxon"
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
B_filteres <- B[,!grepl("Control", colnames(B))]
percent_zero_na <- apply(B_filteres, 1, function(x) mean(x==0|is.na(x))*100)
ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
C=B_filteres+ZZ 
d<-dim(B_filteres)
Cname=colnames(B_filteres)
Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "[.]")
material=rep("NA", Clen)
period=rep("NA",Clen)
sid=rep("NA",Clen)
seqid=rep("NA", Clen)
for(mm in  1:Clen ){
       material[mm]=splitname[[mm]][1]
       period[mm]=splitname[[mm]][2]
       sid[mm]=splitname[[mm]][3]
       seqid[mm]=splitname[[mm]][4]
}
PeriodMaterial=paste0(period, ".", material)
print(paste("genename","genus2","minp"," KP_PeriodMaterial"," KP_period"," KP_material"," KP_sid","  Pwilcox_MC_MP","Pwilcox_MP_OP","Pwilcox_MC_MP","Pwilcox_OC_OP"," Pwilcox_M_O"," Pwilcox_P_C",
            "truefc_MC_MP","truefc_MP_OP","truefc_MC_MP","truefc_OC_OP"," truefc_M_O"," truefc_P_C",
            "meanAll"," meanMC"," meanOC"," meanMP"," meanOP"," meanM","meanC"," meanO"," meanP", "percentZeroNA", sep=","))

for (i in 1:d[1]){
    #print(i)
    genename=A[i,1]
    if(is.na(genename) ){next;} ###remove row with name as NA
    splitG<-strsplit(as.character(genename), "[|.]")
    LLL=length(splitG[[1]])
    genus=splitG[[1]][LLL]
    genus1=str_replace(genus, "[kposgfc]__", "")
    genus2=str_replace(genus1, "_", " ")
    gene<-as.numeric(C[i,])*scaleF
    geneMC=gene[PeriodMaterial=="Menstruating.Cloth"]; meanMC=mean(geneMC)
    geneOC=gene[PeriodMaterial=="Outside_menstruation.Cloth"]; meanOC=mean(geneOC)
    geneMP=gene[PeriodMaterial=="Menstruating.Pad"]; meanMP=mean(geneMP)
    geneOP=gene[PeriodMaterial=="Outside_menstruation.Pad"]; meanOP=mean(geneOP)
    geneM=gene[period=="Menstruating"]; meanM=mean(geneM)
    geneO=gene[period=="Outside_menstruation"]; meanO=mean(geneO)
    geneP=gene[material=="Pad"]; meanP=mean(geneP)
    geneC=gene[material=="Cloth"]; meanC=mean(geneC)

    mydata=data.frame(gene,log(gene),material,period, PeriodMaterial, sid)
    KP_PeriodMaterial=kruskal.test(gene ~ PeriodMaterial, data = mydata)$p.value
    KP_period=kruskal.test(gene ~ period, data = mydata)$p.value
    KP_material=kruskal.test(gene ~ material, data = mydata)$p.value
    KP_sid=kruskal.test(gene ~ sid, data = mydata)$p.value

    meanAll <-mean(gene)
    Pwilcox_MC_OC<-my.wilcox.p.value(geneMC, geneOC, na.rm=TRUE); 
    Pwilcox_MP_OP<-my.wilcox.p.value(geneMP, geneOP, na.rm=TRUE)
    Pwilcox_MC_MP<-my.wilcox.p.value(geneMC, geneMP, na.rm=TRUE)
    Pwilcox_OC_OP<-my.wilcox.p.value(geneOC, geneOP, na.rm=TRUE)
    Pwilcox_M_O<-my.wilcox.p.value(geneM, geneO, na.rm=TRUE)
    Pwilcox_P_C<-my.wilcox.p.value(geneP, geneC, na.rm=TRUE)
    
    truefc_MC_OC=truefc(meanMC/meanOC)
    truefc_MP_OP=truefc(meanMP/meanOP)
    truefc_MC_MP=truefc(meanMC/meanMP)
    truefc_OC_OP=truefc(meanOC/meanOP)
    truefc_M_O=truefc(meanM/meanO)
    truefc_P_C=truefc(meanP/meanC)

   original=B_filteres[i,]
   percentZeroNA=100*(sum(B_filteres[i,]==0)+sum(is.na(B_filteres[i,])))/d[2]
    minp=min(KP_PeriodMaterial, KP_period, KP_material, KP_sid,  Pwilcox_MC_MP,Pwilcox_MP_OP,Pwilcox_MC_MP,Pwilcox_OC_OP, Pwilcox_M_O, Pwilcox_P_C,  na.rm = TRUE)

    print(paste(genename,genus2,minp, KP_PeriodMaterial, KP_period, KP_material, KP_sid,  Pwilcox_MC_MP,Pwilcox_MP_OP,Pwilcox_MC_MP,Pwilcox_OC_OP, Pwilcox_M_O, Pwilcox_P_C,
                truefc_MC_MP,truefc_MP_OP,truefc_MC_MP,truefc_OC_OP, truefc_M_O, truefc_P_C,
                meanAll, meanMC, meanOC, meanMP, meanOP, meanM,meanC, meanO, meanP, percentZeroNA, sep=","))
   # if(minp<=0.05){
    #  p1=ggplot(data=mydata,aes(x=PeriodMaterial, y=gene,color=PeriodMaterial))+ylab(paste(genus2, "RA%")) +geom_boxplot()+ geom_jitter(shape=16, position=position_jitter(0.2))+ theme_classic()+stat_compare_means()+ theme(legend.position = "none")   
    #  p2=ggboxplot(mydata, x = "", y = "gene", color = "asfsGroup", palette = "jco",ylab = paste(genus2, "RA%"), xlab = "asfsGroup")+stat_compare_means()
    # p3=ggplot(mydata, aes(x=asfs, y=gene)) +geom_point(color=asfsGroup)+ geom_smooth(method=lm)+ theme_classic()+ggtitle(paste("Vox1", genus,"Pcor=", sprintf("%.2f",PCor)))
    # ggarrange(p1, p3, ncol = 2, nrow =1 ) %>% ggexport(filename=paste0(genus,".png"), width = 1000, height = 600)
    #}
   
} 
