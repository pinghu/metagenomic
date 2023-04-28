rm(list=ls())
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggpubr)
library(rstatix)

filename="metaphlan.estimated_count.filtered.7"
test <- read.table(filename, header = TRUE, row.names = 1, sep="\t")
test[is.na(test)]<-0
d=dim(test)
X=apply(test, 1, sum)
otu_count=test[(X/d[2] >=10),]
dim(otu_count)
test_d=data.frame(t(otu_count))
dim(test_d)
min(apply(test_d, 1, sum)) ###301
test_a <- decostand(test_d, method = "total")
write.csv(t(test_a), file =paste0(filename, "recal.Percent.csv"))

Cname=rownames(test_d)

Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "[.]")


Treat=rep("NA", Clen)
sid=rep("NA", Clen)
Visit=rep("NA", Clen)
TreatVisit=rep("NA", Clen)
asfs=rep("NA", Clen)
age=rep("NA", Clen)
id2=rep("NA", Clen)
asfsGroup=rep("NA",Clen)
ageGroup=rep("NA", Clen)
for(mm in  1:Clen ){
  Treat[mm]=splitname[[mm]][2]
  Visit[mm]=splitname[[mm]][3]
  sid[mm]=splitname[[mm]][4]
  asfs[mm]=as.numeric(splitname[[mm]][5])
  age[mm]=splitname[[mm]][6]
  id2[mm]=splitname[[mm]][7]
  if(age[mm]<40){ageGroup[mm]="40Less"}
  else if (age[mm]<50){ageGroup[mm]="50-40"}
  else if (age[mm]<60){ageGroup[mm]="60-50"}
  else{ageGroup[mm]="70+-60"}
  if(asfs[mm]>=24){asfsGroup[mm]="3D"}
  else if(asfs[mm]<=10){asfsGroup[mm]="1ND"}
  else{asfsGroup[mm]="2MID"}
}
Visit[Visit=="BL"]="Baseline"
Visit[Visit=="W3"]="Week3"
Treat[Treat=="B"]="Piroctone_Olamine"
Treat[Treat=="C"]="Control"
Treat[Treat=="A"]="ZPT"
TreatVisit=paste0(Treat, Visit)
shannon<-diversity(test_d)
simp<-diversity(test_d, "simpson")
invsimp<-diversity(test_d, "inv")
#unbias.simp <- rarefy(test_d, 2) - 1
observed<-apply(test_d>0,1,sum)
N <- apply(test_d,1,sum)
Menhinick_index<-observed/sqrt(N)
Margalef_index <-(observed-1)/log(N)
mydata=data.frame(shannon,simp,invsimp,  observed, Menhinick_index, Margalef_index, 
                  as.numeric(asfs), asfsGroup, TreatVisit,Treat, Visit, age, ageGroup,sid, id2 )

write.csv(t(mydata), file =paste0(filename, "alphadiversity.csv"))




png(filename=paste0(filename,"shannon.TreatVisit.wilcox.png"), width = 700, height =600)
stat.test <- mydata %>%
  #group_by(Country) %>%
  wilcox_test(shannon ~ TreatVisit) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 2)
bxp <- ggboxplot(mydata, x = "TreatVisit", y = "shannon", fill = "Treat")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+ theme(legend.position = "none")+geom_jitter()
dev.off()




png(filename=paste0(filename,"observed.TreatVisit.wilcox.png"), width = 700, height = 600)
stat.test <- mydata %>%
  wilcox_test(observed ~ TreatVisit) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 15.8)
bxp <- ggboxplot(mydata, x = "TreatVisit", y = "observed", fill = "Treat")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+ theme(legend.position = "none")+geom_jitter()
dev.off()

KP_asfsGroup=kruskal.test(shannon ~ asfsGroup, data = mydata)$p.value
KP_TreatVisit=kruskal.test(shannon ~TreatVisit, data = mydata)$p.value
KP_ageGroup=kruskal.test(shannon ~ ageGroup, data = mydata)$p.value
KP_Treat=kruskal.test(shannon ~Treat, data = mydata)$p.value
KP_Visit=kruskal.test(shannon ~Visit, data = mydata)$p.value
print(paste0("Shannon Diversity: KP_asfsGroup=", KP_asfsGroup, "; KP_TreatVisit=", KP_TreatVisit,"; KP_Treat=", KP_Treat,"; KP_Visit=", KP_Visit, "; KP_ageGroup=", KP_ageGroup
))


KP_asfsGroup2=kruskal.test(observed ~ asfsGroup, data = mydata)$p.value
KP_TreatVisit2=kruskal.test(observed ~TreatVisit, data = mydata)$p.value
KP_ageGroup2=kruskal.test(observed ~ ageGroup, data = mydata)$p.value
KP_Treat2=kruskal.test(observed ~Treat, data = mydata)$p.value
KP_Visit2=kruskal.test(observed ~Visit, data = mydata)$p.value
print(paste0("Oberved: KP_asfsGroup=", KP_asfsGroup2, "; KP_TreatVisit=", KP_TreatVisit2,"; KP_Treat=", KP_Treat2,"; KP_Visit=", KP_Visit2, "; KP_ageGroup=", KP_ageGroup2
))
#[1] "Shannon Diversity: KP_asfsGroup=0.00636732064441896; KP_TreatVisit=7.81609751196999e-08; KP_Treat=0.000208022003009667; KP_Visit=0.00375713990853845; KP_ageGroup=0.0164831680750766"
#[1] "Oberved: KP_asfsGroup=0.846354807337587; KP_TreatVisit=0.638978033915691; KP_Treat=0.419638290070914; KP_Visit=0.359382764872769; KP_ageGroup=2.51569355714698e-08"

x<-test_d
beta_dist <- vegdist(x, index = "bray")
mds <- metaMDS(beta_dist)
mds_data <- as.data.frame(mds$points)
mds_data$SampleID <- rownames(mds_data)
splitname<-strsplit(mds_data$SampleID, "[._]")


mds_data$asfs=rep("NA", Clen)
mds_data$asfsGroup=rep("NA", Clen)
mds_data$age=rep("NA", Clen)
mds_data$ageGroup=rep("NA", Clen)
mds_data$sid=rep("NA", Clen)
mds_data$TreatVisit=rep("NA", Clen)
mds_data$Visit=rep("NA", Clen)
mds_data$Treat=rep("NA", Clen)
mds_data$id2=rep("NA", Clen)
splitname<-strsplit(Cname, "[._]")
for(mm in  1:Clen ){
  mds_data$Treat[mm]=splitname[[mm]][2]
  mds_data$Visit[mm]=splitname[[mm]][3]
  mds_data$sid[mm]=splitname[[mm]][4]
  mds_data$asfs[mm]=as.numeric(splitname[[mm]][5])
  mds_data$age[mm]=splitname[[mm]][6]
  mds_data$id2[mm]=splitname[[mm]][7]
  
  if(mds_data$age[mm]<40){mds_data$ageGroup[mm]="40Less"}
  else if (mds_data$age[mm]<50){mds_data$ageGroup[mm]="50-40"}
  else if (mds_data$age[mm]<60){mds_data$ageGroup[mm]="60-50"}
  else{mds_data$ageGroup[mm]="70+-60"}
  if(mds_data$asfs[mm]>=24){mds_data$asfsGroup[mm]="3D"}
  else if(mds_data$asfs[mm]<=10){mds_data$asfsGroup[mm]="1ND"}
  else{mds_data$asfsGroup[mm]="2MID"}
  
}
mds_data$Visit[mds_data$Visit=="BL"]="Baseline"
mds_data$Visit[mds_data$Visit=="W3"]="Week3"
mds_data$Treat[mds_data$Treat=="B"]="Piroctone_Olamine"
mds_data$Treat[mds_data$Treat=="C"]="Control"
mds_data$Treat[mds_data$Treat=="A"]="ZPT"
mds_data$TreatVisit=paste0(mds_data$Treat, mds_data$Visit)

bc.dist <- vegdist(test_d, method = "bray")
png(filename=paste0(filename,".beta.Visit.facetTreat.png"),  width=600, height=380)
p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Visit)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=Visit, fill=Visit), alpha=0.1, geom="polygon")+ facet_grid(. ~ Treat)
print(p1)
dev.off()

print("adonis result bc.dist ~ TreatVisit")
adonis2(bc.dist ~ TreatVisit)
test.adonis <- adonis(bc.dist ~ TreatVisit)
test.adonis <- as.data.frame(test.adonis$aov.tab)
test.adonis
#Df SumsOfSqs   MeanSqs  F.Model        R2 Pr(>F)
#TreatVisit   3  3.854484 1.2848281 7.409315 0.1018565  0.001
#Residuals  196 33.987797 0.1734071       NA 0.8981435     NA
#Total      199 37.842281        NA       NA 1.0000000     NA
print("Pairwise Anosim Test")
cbn <- combn(x=unique(TreatVisit), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- test_d [c(TreatVisit %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(TreatVisit %in% cbn[,i]), ]
  
  permanova_pairwise <- adonis(vegdist(ps.subs, method = "bray") ~ metadata_sub$TreatVisit)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
p.table

#1                      2                           p     p.adj
#1 Piroctone_OlamineBaseline Piroctone_OlamineWeek3 0.001 0.0020
#2 Piroctone_OlamineBaseline        ControlBaseline 0.006 0.0090
#3 Piroctone_OlamineBaseline           ControlWeek3 0.123 0.1476
#4    Piroctone_OlamineWeek3        ControlBaseline 0.001 0.0020
#5    Piroctone_OlamineWeek3           ControlWeek3 0.001 0.0020
#6           ControlBaseline           ControlWeek3 0.236 0.2360

