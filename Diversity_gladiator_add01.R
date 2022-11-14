rm(list=ls())
#library(tidyverse)
# p1 <- c("tidyverse", "vegan", "BiocManager")
# p2 <- c("phyloseq", "ANCOMBC", "DESeq2", "ComplexHeatmap")
# load_package <- function(p) {
#   if (!requireNamespace(p, quietly = TRUE)) {
#     ifelse(p %in% p1, 
#            install.packages(p, repos = "http://cran.us.r-project.org/"), 
#            BiocManager::install(p))
#   }
#   library(p, character.only = TRUE, quietly = TRUE)
# }
# invisible(lapply(c(p1,p2), load_package))

args <- commandArgs(trailingOnly = TRUE)
filename<- args[1]
#filename="gladiator.metaphlan3.count10.7"
#filename="gladiator.metaphlan3.count10.7.bacteria"
filename="gladiator.metaphlan3.count10.7.fungal"  
#filename="gladiator.kraken.count.10.species"
#-- could not do beta analysis as some empty can I add 1 to all 

test <- read.table(filename, header = TRUE, row.names = 1, sep="\t")
test[is.na(test)]<-0
d=dim(test)
d
X=apply(test, 1, sum)
test_filter2=test[(X/d[2] >=10),]

test_d=data.frame(t(test_filter2))+0.01
dim(test_d)
min(apply(test_d, 1, sum)) ###301


library(vegan)
# Turn percent cover to relative abundance by dividing each value by sample
# total abundance
test_a <- decostand(test_d, method = "total")
# check total abundance in each sample
#apply(test_a, 1, sum) ### this is just for confirmation
#------------------------------------------
write.table(t(test_a), file =paste0(filename, "recal.Percent.RShort"), col.names=TRUE, row.names=TRUE, sep = "\t")

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
  TreatVisit[mm]=paste0(Treat[mm], Visit[mm])
  if(age[mm]<40){ageGroup[mm]="40Less"}
  else if (age[mm]<50){ageGroup[mm]="50-40"}
  else if (age[mm]<60){ageGroup[mm]="60-50"}
  else{ageGroup[mm]="70+-60"}
  if(asfs[mm]>=24){asfsGroup[mm]="3D"}
  else if(asfs[mm]<=10){asfsGroup[mm]="1ND"}
  else{asfsGroup[mm]="2MID"}
}



library(stringr)

library(ggplot2)

#https://stackoverflow.com/questions/30057765/histogram-ggplot-show-count-label-for-each-bin-for-each-category

shannon<-diversity(test_d)
simp<-diversity(test_d, "simpson")
invsimp<-diversity(test_d, "inv")
## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
#unbias.simp <- rarefy(test_d, 2) - 1
## Fisher alpha
#fisher <- fisher.alpha(test_d)
####observed species
observed<-apply(test_d>0,1,sum)
N <- apply(test_d,1,sum)
###Richness Index#####
Menhinick_index<-observed/sqrt(N)
Margalef_index <-(observed-1)/log(N)
#unbias.simp
mydata=data.frame(shannon,simp,invsimp,  observed, Menhinick_index, Margalef_index, 
                   as.numeric(asfs), asfsGroup, TreatVisit,Treat, Visit, age, ageGroup,sid, id2 )

write.table(t(mydata), file =paste0(filename, "alphadiversity.RShort"), col.names=TRUE, row.names=TRUE, sep = "\t")

saveRDS(mydata,paste0(filename, ".diversitydata"))
saveRDS(test_d,paste0(filename,".test_d"))



library(ggplot2)

#https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/
library(tidyverse)
library(ggpubr)
library(rstatix)

png(filename=paste0(filename,"shannon.asfsGroup.ttest.png"))
stat.test <- mydata %>%
  #group_by(Country) %>%
  t_test(shannon ~ asfsGroup) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 3.8)
bxp <- ggboxplot(mydata, x = "asfsGroup", y = "shannon", color = "asfsGroup")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
dev.off()

png(filename=paste0(filename,"shannon.asfsGroup.wilcox.png"))
stat.test <- mydata %>%
  #group_by(Country) %>%
  wilcox_test(shannon ~ asfsGroup) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 3.8)
bxp <- ggboxplot(mydata, x = "asfsGroup", y = "shannon", color = "asfsGroup")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
dev.off()

png(filename=paste0(filename,"shannon.TreatVisit.wilcox.png"))
stat.test <- mydata %>%
  #group_by(Country) %>%
  wilcox_test(shannon ~ TreatVisit) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 3.8)
bxp <- ggboxplot(mydata, x = "TreatVisit", y = "shannon", color = "TreatVisit")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
dev.off()

wilcox.shannon <- pairwise.wilcox.test(mydata$shannon, 
                                        mydata$TreatVisit, 
                                        p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("shannon -- TreatVisit wilcoxon adjustp")
tab.shannon

png(filename=paste0(filename,"shannon.TreatVisit.ttest.png"))
stat.test <- mydata %>%
  #group_by(Country) %>%
  t_test(shannon ~ TreatVisit) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 3.8)
bxp <- ggboxplot(mydata, x = "TreatVisit", y = "shannon", color = "TreatVisit")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
dev.off()

# png(filename=paste0(filename,"observed.TreatVisit.ttest.png"))
# stat.test <- mydata %>%
#   t_test(observed ~ TreatVisit) %>%
#   #adjust_pvalue() %>%
#   mutate(y.position = 15.8)
# bxp <- ggboxplot(mydata, x = "TreatVisit", y = "observed", color = "TreatVisit")
# bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
# dev.off()

# png(filename=paste0(filename,"observed.TreatVisit.wilcox.png"))
# stat.test <- mydata %>%
#   wilcox_test(observed ~ TreatVisit) %>%
#   #adjust_pvalue() %>%
#   mutate(y.position = 15.8)
# bxp <- ggboxplot(mydata, x = "TreatVisit", y = "observed", color = "TreatVisit")
# bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
# dev.off()
########add in to show statistical table directly
# wilcox.observed <- pairwise.wilcox.test(mydata$observed, 
#                                         mydata$TreatVisit, 
#                                         p.adjust.method = "BH")
# tab.observed <- wilcox.observed$p.value %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "group1") %>%
#   gather(key="group2", value="p.adj", -group1) %>%
#   na.omit()
# print("observed -- TreatVisit wilcoxon adjustp")
# tab.observed


# png(filename=paste0(filename,"shannon.Visit.png"))
# stat.test <- mydata %>%
#   #group_by(Country) %>%
#   t_test(shannon ~ Visit) %>%
#   #adjust_pvalue() %>%
#   mutate(y.position = 3.8)
# bxp <- ggboxplot(mydata, x = "Visit", y = "shannon", color = "Visit")
# bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
# dev.off()


# png(filename=paste0(filename,"shannon.Treat.png"))
# stat.test <- mydata %>%
#   #group_by(Country) %>%
#   t_test(shannon ~ Treat) %>%
#   #adjust_pvalue() %>%
#   mutate(y.position = 3.8)
# bxp <- ggboxplot(mydata, x = "Treat", y = "shannon", color = "Treat")
# bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
# dev.off()

# 
# 
# png(filename=paste0(filename,"shannon.ageGroup.png"))
# stat.test <- mydata %>%
#   #group_by(Country) %>%
#   t_test(shannon ~ ageGroup) %>%
#   #adjust_pvalue() %>%
#   mutate(y.position = 3.8)
# bxp <- ggboxplot(mydata, x = "ageGroup", y = "shannon", color = "ageGroup")
# bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
# dev.off()
colnames(mydata)[7]="asfs"
png(filename=paste0(filename,"asfs.TreatVisit.ttest.png"))
stat.test <- mydata %>%
  #group_by(Country) %>%
  t_test(asfs ~ TreatVisit) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 55)
bxp <- ggboxplot(mydata, x = "TreatVisit", y = "asfs", color = "TreatVisit")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
dev.off()

png(filename=paste0(filename,"asfs.TreatVisit.wilcox.png"))
stat.test <- mydata %>%
  #group_by(Country) %>%
  wilcox_test(asfs ~ TreatVisit) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 55)
bxp <- ggboxplot(mydata, x = "TreatVisit", y = "asfs", color = "TreatVisit")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
dev.off()
wilcox.asfs <- pairwise.wilcox.test(mydata$asfs, 
                                        mydata$TreatVisit, 
                                        p.adjust.method = "BH")
tab.asfs <- wilcox.asfs$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("asfs -- TreatVisit wilcoxon adjustp")
tab.asfs



KP_asfsGroup=kruskal.test(shannon ~ asfsGroup, data = mydata)$p.value
KP_TreatVisit=kruskal.test(shannon ~TreatVisit, data = mydata)$p.value
KP_ageGroup=kruskal.test(shannon ~ ageGroup, data = mydata)$p.value
KP_Treat=kruskal.test(shannon ~Treat, data = mydata)$p.value
KP_Visit=kruskal.test(shannon ~Visit, data = mydata)$p.value
print(paste0("KP_asfsGroup=", KP_asfsGroup, "; KP_TreatVisit=", KP_TreatVisit,"; KP_Treat=", KP_Treat,"; KP_Visit=", KP_Visit, "; KP_ageGroup=", KP_ageGroup
     ))
####bacteria and fungal [1] "KP_asfsGroup=5.72e-10; KP_TreatVisit=9.74e-27; KP_Treat=3.78e-13; KP_Visit=1.23e-11; KP_ageGroup=0.01"
#####bacteria and fungal +0.01 : [1] "KP_asfsGroup=4.99e-10; KP_TreatVisit=9.44e-27; KP_Treat=3.65e-13; KP_Visit=1.24e-11; KP_ageGroup=0.006"
#####https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html#######There is code for all the diversity calculation
#####Bacteria_only result:[1] "KP_asfsGroup=0.45; KP_TreatVisit=0.20; KP_Treat=0.29; KP_Visit=0.14; KP_ageGroup=4.35e-07"
#####Bacteria +0.01 [1] "KP_asfsGroup=0.30; KP_TreatVisit=0.08; KP_Treat=0.16; KP_Visit=0.10; KP_ageGroup=2.65e-07"
#####Bacteria+1: [1] "KP_asfsGroup=0.25; KP_TreatVisit=0.053; KP_Treat=0.09; KP_Visit=0.12; KP_ageGroup=2.88e-07"
#####Bacteria+0.1 [1] "KP_asfsGroup=0.30; KP_TreatVisit=0.08; KP_Treat=0.14; KP_Visit=0.12; KP_ageGroup=2.61e-07"
#####Fungal only result: [1] "KP_asfsGroup=9.04e-08; KP_TreatVisit=1.41e-22; KP_Treat=2.34e-10; KP_Visit=7.52e-12; KP_ageGroup=0.60"
#####Fungal add 0.01 [1] "KP_asfsGroup=0.002; KP_TreatVisit=0.03; KP_Treat=0.43; KP_Visit=0.03; KP_ageGroup=0.35"
#####Fungal add 0.1 [1] "KP_asfsGroup=0.002; KP_TreatVisit=0.03; KP_Treat=0.43; KP_Visit=0.03; KP_ageGroup=0.36"
#####Fungal Only add 1 to matrix: [1] "KP_asfsGroup=0.002; KP_TreatVisit=0.03; KP_Treat=0.43; KP_Visit=0.03; KP_ageGroup=0.36"
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
  mds_data$TreatVisit[mm]=paste0(mds_data$Treat[mm], mds_data$Visit[mm])
  if(mds_data$age[mm]<40){mds_data$ageGroup[mm]="40Less"}
  else if (mds_data$age[mm]<50){mds_data$ageGroup[mm]="50-40"}
  else if (mds_data$age[mm]<60){mds_data$ageGroup[mm]="60-50"}
  else{mds_data$ageGroup[mm]="70+-60"}
  if(mds_data$asfs[mm]>=24){mds_data$asfsGroup[mm]="3D"}
  else if(mds_data$asfs[mm]<=10){mds_data$asfsGroup[mm]="1ND"}
  else{mds_data$asfsGroup[mm]="2MID"}
  
}

library(stringr)

#install.packages('devtools')
#library(devtools)
#install_github('fawda123/ggord')
#library(ggord)
library(ggplot2)
 
 png(filename=paste0(filename,".beta.asfsGroup.png"),  width=600, height=380)
# p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = TreatCode)) +geom_point()+theme_bw()+stat_ellipse()
 #+geom_text(aes(label=SID),hjust=0.5, vjust=0)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = asfsGroup)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=asfsGroup, fill=asfsGroup), alpha=0.1, geom="polygon")
 print(p1)
 dev.off()
 
 png(filename=paste0(filename,".beta.asfsGroup.facetTreat.png"),  width=600, height=380)
 # p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = TreatCode)) +geom_point()+theme_bw()+stat_ellipse()
 #+geom_text(aes(label=SID),hjust=0.5, vjust=0)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = asfsGroup)) +geom_point(size=3)+theme_bw()+stat_chull(aes(color=asfsGroup, fill=asfsGroup), alpha=0.1, geom="polygon")+ facet_grid(. ~ Treat)
 print(p1)
 dev.off()
 
 # png(filename=paste0(filename,".beta.TreatVisit-facet.png"),  width=600, height=380)
 # p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = TreatVisit)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=TreatVisit, fill=TreatVisit), alpha=0.1, geom="polygon")+ facet_grid(. ~ Treat)
 # print(p1)
 # dev.off()
 png(filename=paste0(filename,".beta.Visit.facetTreat.png"),  width=600, height=380)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Visit)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=Visit, fill=Visit), alpha=0.1, geom="polygon")+ facet_grid(. ~ Treat)
 print(p1)
 dev.off()
 # png(filename=paste0(filename,".beta.Visit.png"),  width=600, height=380)
 # p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Visit)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=Visit, fill=Visit), alpha=0.1, geom="polygon")
 # print(p1)
 # dev.off()
 # png(filename=paste0(filename,".beta.Treat.png"),  width=600, height=380)
 # p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Treat)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=Treat, fill=Treat), alpha=0.1, geom="polygon")
 # print(p1)
 # dev.off()
 # 
 # 
 # png(filename=paste0(filename,".beta.ageGroup.png"),  width=600, height=380)
 # #p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Visit)) +geom_point()+theme_bw()+stat_ellipse()
 # #+geom_text(aes(label=SID),hjust=0.5, vjust=0)
 # p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = ageGroup)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=ageGroup, fill=ageGroup), alpha=0.1, geom="polygon")
 # 
 # print(p1)
 # dev.off()



############################################################
print("all samples")
# calculate Bray-Curtis distance among samples

# ps.rarefied=rarify(test_d)
# dist = phyloseq::distance(ps.rarefied, method="bray")
# ordination = ordinate(ps.rarefied, method="PCoA", distance=dist)
# plot_ordination(ps.rarefied, ordination, color="TreatVisit") + 
#   theme_classic() +
#   theme(strip.background = element_blank())
# metadata <- data.frame(sample_data(ps.rarefied))
# test.adonis <- adonis(dist ~ body.site, data = metadata)
# test.adonis <- as.data.frame(test.adonis$aov.tab)
# test.adonis
 ######################################
 
 ###################################################################### 
# anosim(dist, metadata$body.site) 
#####what is in mds_data
bc.dist <- vegdist(test_d, method = "bray")
jaccard.dist <- vegdist(test_d, method = "jaccard")
mds2 <- metaMDS(jaccard.dist)
mds_data2 <- as.data.frame(mds2$points)
mds_data2$SampleID <- rownames(mds_data2)
mds_data2$asfs=mds_data$asfs
mds_data2$asfsGroup=mds_data$asfsGroup
mds_data2$age=mds_data$age
mds_data2$ageGroup=mds_data$ageGroup
mds_data2$sid=mds_data$sid
mds_data2$TreatVisit=mds_data$TreatVisit
mds_data2$Visit=mds_data$Visit
mds_data2$Treat=mds_data$Treat
mds_data2$id2=mds_data$id2
png(filename=paste0(filename,".beta.Visit.facetTreat.jaccard.png"),  width=600, height=380)
p1<-ggplot(mds_data2, aes(x = MDS1, y = MDS2, color = Visit)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=Visit, fill=Visit), alpha=0.1, geom="polygon")+ facet_grid(. ~ Treat)
print(p1)
dev.off()
# cluster communities using average-linkage algorithm
bc.clust <- hclust(bc.dist, method = "average")
# plot cluster diagram
png(filename=paste0(filename,".bray-Curtis-cluster.png"), width=8000, height=2000)
plot(bc.clust, ylab = "Bray-Curtis dissimilarity")
dev.off()

print("adonis result bc.dist ~ TreatVisit")
adonis2(bc.dist ~ TreatVisit)
#bacteriaonly p=0.001
#bacteria+0.01 p=0.001
#bacteria and fungi p=0.001
#fungal+1 p=0.001
test.adonis <- adonis(bc.dist ~ TreatVisit)
test.adonis <- as.data.frame(test.adonis$aov.tab)
test.adonis

#PAIRWISE PERMANOVA
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
###########################Jaccard test
print("adonis result jaccard.dist ~ TreatVisit")
adonis2(jaccard.dist ~ TreatVisit)
#bacteriaonly p=0.001
#bacteria+0.01 p=0.001
#bacteria and fungi p=0.001
#fungal+1 p=0.001
test.adonis <- adonis(jaccard.dist ~ TreatVisit)
test.adonis <- as.data.frame(test.adonis$aov.tab)
test.adonis

#PAIRWISE PERMANOVA
cbn <- combn(x=unique(TreatVisit), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- test_d [c(TreatVisit %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(TreatVisit %in% cbn[,i]), ]
  
  permanova_pairwise <- adonis(vegdist(ps.subs, method = "jaccard") ~ metadata_sub$TreatVisit)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
p.table

#####################################################
####anosim test
print("Anosim Test of Jaccard distance to Treatvisit")
anosim(jaccard.dist, TreatVisit)
print("Pairwise Anosim Test")
cbn <- combn(x=unique(TreatVisit), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- test_d [c(TreatVisit %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(TreatVisit %in% cbn[,i]), ]
  
  permanova_pairwise <- anosim(vegdist(ps.subs, method = "jaccard") , metadata_sub$TreatVisit)
  p <- c(p, permanova_pairwise$signif[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
p.table



#print("adonis result bc.dist ~ Treat")
#adonis2(bc.dist ~ Treat)
#bacteria only p=0.001
# bacteria and fungi p=0.001
#fungal+1 p=0.001

#print("adonis result bc.dist ~ Visit")
#adonis2(bc.dist ~ Visit)
#bacteria only p=0.001
#bacteri and fungi p=0.001
#fungal+1 p=0.001
# Taxonomic (Bray-Curtis) dissimilarity explained
#print("adonis result bc.dist ~ asfsGroup")
#adonis2(bc.dist ~ asfsGroup)


#bacteriaonly p=0.001
#bacteriaonly+0.01 p=0.001
#bacteria+fungal p=0.001
#fungal+1 p=0.001
#print("adonis  result bc.dist ~ ageGroup")
#adonis2(bc.dist ~ ageGroup)
# bacteria only p=0.175
#bacteria + fungal p=0.173
#fungal+1 p=0.336


  
############################################################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("M3C")
#library(M3C)
##########################################################################
#https://www.molecularecologist.com/2013/08/20/making-heatmaps-with-r-for-microbiome-analysis/
#https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01013-0 
###############################################################################################
###Now look into the bar plot and the enrichment analysis
###########################################
