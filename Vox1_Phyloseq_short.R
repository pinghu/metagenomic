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

asfs=rep("NA", Clen)
asfsGroup=rep("NA",Clen)
gender=rep("NA",Clen)
sid=rep("NA",Clen)
age=rep("NA", Clen)
ageGroup=rep("NA", Clen)
category=rep("NA", Clen)
for(mm in  1:Clen ){
  asfs[mm]=as.numeric(splitname[[mm]][3])
  category[mm]=splitname[[mm]][2]
  age[mm]=splitname[[mm]][5]
  gender[mm]=splitname[[mm]][4]
  sid[mm]=splitname[[mm]][6]
  if(age[mm]<40){ageGroup[mm]="40Less"}
  else if (age[mm]<50){ageGroup[mm]="50-40"}
  else if (age[mm]<60){ageGroup[mm]="60-50"}
  else{ageGroup[mm]="70+-60"}
}
asfs=as.numeric(asfs)
asfsGroup[asfs>=24]="3Dandruff"
asfsGroup[asfs<10]="1NoDandruff"
asfsGroup[asfs<24&asfs>=10]="2Intermediate"
category == asfsGroup

age=as.numeric(age)
ageCategory2=rep("NA", Clen)
ageCategory2[age<40]="40Less"
ageCategory2[age<50 & age>=40]="50-40"
ageCategory2[age<60 & age >=50]="60-50"
ageCategory2[age>=60]="70+-60"

ageGroup ==ageCategory2


shannon<-diversity(test_d)
simp<-diversity(test_d, "simpson")
invsimp<-diversity(test_d, "inv")
#unbias.simp <- rarefy(test_d, 2) - 1
observed<-apply(test_d>0,1,sum)
N <- apply(test_d,1,sum)
Menhinick_index<-observed/sqrt(N)
Margalef_index <-(observed-1)/log(N)
mydata=data.frame(shannon,simp,invsimp,  observed, Menhinick_index, Margalef_index, 
                  asfs, asfsGroup, category, age, ageGroup,sid, gender )

write.csv(t(mydata), file =paste0(filename, "alphadiversity.csv"))


png(filename=paste0(filename,"shannon.asfsGroup.wilcox.png"), width=300, height = 600)
stat.test <- mydata %>%
  #group_by(Country) %>%
  wilcox_test(shannon ~ asfsGroup) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 2)
bxp <- ggboxplot(mydata, x = "asfsGroup", y = "shannon", fill = "asfsGroup")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+ theme(legend.position = "none")+geom_jitter()
dev.off()


png(filename=paste0(filename,"observed.Category.wilcox.png"), width=300, height = 380)
stat.test <- mydata %>%
  wilcox_test(observed ~ asfsGroup) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 40)
bxp <- ggboxplot(mydata, x = "asfsGroup", y = "observed", fill = "asfsGroup")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+ theme(legend.position = "none")+geom_jitter()
dev.off()

KP_asfsGroup=kruskal.test(shannon ~ asfsGroup, data = mydata)$p.value
KP_Category=kruskal.test(shannon ~category, data = mydata)$p.value
KP_ageGroup=kruskal.test(shannon ~ ageGroup, data = mydata)$p.value
print(paste0("Shannon KP_asfsGroup=", KP_asfsGroup, "; KP_Category=", KP_Category, "; KP_ageGroup=", KP_ageGroup
))
KP_asfsGroup2=kruskal.test(observed ~ asfsGroup, data = mydata)$p.value
KP_Category2=kruskal.test(observed ~category, data = mydata)$p.value
KP_ageGroup2=kruskal.test(observed ~ ageGroup, data = mydata)$p.value
print(paste0("observed KP_asfsGroup2=", KP_asfsGroup2, "; KP_Category2=", KP_Category2, "; KP_ageGroup2=", KP_ageGroup2
))


#[1] "observed KP_asfsGroup2=0.185; KP_Category2=0.185; KP_ageGroup2=0.339"
#[1] "Shannon KP_asfsGroup=0.708; KP_Category=0.708; KP_ageGroup=0.356"

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
mds_data$gender=rep("NA", Clen)
splitname<-strsplit(Cname, "[._]")
for(mm in  1:Clen ){
  mds_data$asfs[mm]=as.numeric(splitname[[mm]][3])
  mds_data$asfsGroup[mm]=splitname[[mm]][2]
  mds_data$gender[mm]=splitname[[mm]][4]
  mds_data$age[mm]=splitname[[mm]][5]
  mds_data$sid[mm]=splitname[[mm]][6]
  
  if(mds_data$age[mm]<40){mds_data$ageGroup[mm]="40Less"}
  else if (mds_data$age[mm]<50){mds_data$ageGroup[mm]="50-40"}
  else if (mds_data$age[mm]<60){mds_data$ageGroup[mm]="60-50"}
  else{mds_data$ageGroup[mm]="70+-60"}
  
}
mds_data$asfs=as.numeric(mds_data$asfs)
mds_data$age=as.numeric(mds_data$age)
mds_data$asfsGroup[mds_data$asfsGroup=="1ND"]="1NoDandruff"
mds_data$asfsGroup[mds_data$asfsGroup=="2MID"]="2Intermediate"
mds_data$asfsGroup[mds_data$asfsGroup=="3D"]="3Dandruff"
png(filename=paste0(filename,".beta.asfsGroup.png"),  width=600, height=380)
p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = asfsGroup)) +geom_point(size=5)+theme_bw()+stat_chull(aes(color=asfsGroup, fill=asfsGroup), alpha=0.1, geom="polygon")
print(p1)
dev.off()

#########################################################
bc.dist <- vegdist(test_d, method = "bray")


print("adonis result bc.dist ~ asfsGroup")
adonis2(bc.dist ~ asfsGroup)
test.adonis <- adonis(bc.dist ~ asfsGroup)
test.adonis <- as.data.frame(test.adonis$aov.tab)
test.adonis
#adonis2(formula = bc.dist ~ asfsGroup)
#Df SumOfSqs      R2      F Pr(>F)    
#asfsGroup  2   2.0803 0.09851 4.9722  0.001 ***
#  Residual  91  19.0366 0.90149                  
#Total     93  21.1169 1.00000  

########################################
cbn <- combn(x=unique(asfsGroup), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- test_d [c(asfsGroup %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(asfsGroup %in% cbn[,i]), ]
  
  permanova_pairwise <- adonis(vegdist(ps.subs, method = "bray") ~ metadata_sub$asfsGroup)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
p.table
#1             2     p  p.adj
#1   1NoDandruff 2Intermediate 0.566 0.5660
#2   1NoDandruff     3Dandruff 0.001 0.0015
#3 2Intermediate     3Dandruff 0.001 0.0015
