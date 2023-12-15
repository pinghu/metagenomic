rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
filename<- args[1]
filename="MulSanne.species10.txt"

#-- could not do beta analysis as some empty can I add 1 to all 

test <- read.table(filename, header = TRUE, row.names = 1, sep="\t")
test[is.na(test)]<-0
d=dim(test)
d
X=apply(test, 1, sum)
test_filter2=test[(X/d[2] >=10),]

test_d=data.frame(t(test_filter2))
dim(test_d)
min(apply(test_d, 1, sum))


library(vegan)
# Turn percent cover to relative abundance by dividing each value by sample
# total abundance
test_a <- decostand(test_d, method = "total")
# check total abundance in each sample
#apply(test_a, 1, sum) ### this is just for confirmation
#------------------------------------------
write.table(t(test_a), file =paste0(filename, ".recal.Percent.RShort"), col.names=TRUE, row.names=TRUE, sep = "\t")

Cname=rownames(test_d)

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
  "HEALTHY_BL", "UNHEALTHY_BL", "HEALTHY_W4","UNHEALTHY_W4",  "HEALTHY_W8",   "UNHEALTHY_W8"    
  
)
DiseaseVisit <- factor(DiseaseVisit, levels = custom_order)



library(stringr)

library(ggplot2)

#https://stackoverflow.com/questions/30057765/histogram-ggplot-show-count-label-for-each-bin-for-each-category
shannon<-diversity(test_d)
observed<-diversity(test_d)
simp<-diversity(test_d, "simpson")
invsimp<-diversity(test_d, "inv")
## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
integer_df <- apply(test_d, 2, as.integer)
unbias.simp <- rarefy(integer_df, 2) - 1
## Fisher alpha
#fisher <- fisher.alpha(test_d)
####observed species
observed<-apply(test_d>0,1,sum)
N <- apply(test_d,1,sum)
###Richness Index#####
Menhinick_index<-observed/sqrt(N)
Margalef_index <-(observed-1)/log(N)
#unbias.simp
mydata=data.frame(shannon, observed,simp,invsimp,  observed, Menhinick_index, Margalef_index, Visit, Disease, SID, Gender, DiseaseVisit, DiseaseVisitSID)
                

write.table(t(mydata), file =paste0(filename, ".alphadiversity.RShort"), col.names=TRUE, row.names=TRUE, sep = "\t")

saveRDS(mydata,paste0(filename, ".diversitydata"))
saveRDS(test_d,paste0(filename,".test_d"))



library(ggplot2)

#https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/
library(tidyverse)
library(ggpubr)
library(rstatix)

png(filename=paste0(filename,".observed.DiseaseVisit.ttest.png"), width=800, height = 800)
stat.test <- mydata %>%
  #group_by(Country) %>%
  t_test(observed ~ DiseaseVisit) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 160)
bxp <- ggboxplot(mydata, x = "DiseaseVisit", y = "observed", color = "Disease")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
dev.off()


png(filename=paste0(filename,".observed.DiseaseVisit.wilcox.png"), width=800, height = 800)
stat.test <- mydata %>%
  #group_by(Country) %>%
  wilcox_test(observed ~ DiseaseVisit) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 160)
bxp <- ggboxplot(mydata, x = "DiseaseVisit", y = "observed", color = "Disease")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
dev.off()

wilcox.observed <- pairwise.wilcox.test(mydata$observed, 
                                        mydata$DiseaseVisit, 
                                        p.adjust.method = "BH")
tab.observed <- wilcox.observed$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("observed -- DiseaseVisit wilcoxon adjustp")
tab.observed

write.table(tab.observed, file = "MulSanne.observed.wilcox.p.txt", sep = "\t", row.names = FALSE)
#############################


png(filename=paste0(filename,".shannon.DiseaseVisit.ttest.png"), width=800, height = 800)
stat.test <- mydata %>%
  #group_by(Country) %>%
  t_test( shannon~ DiseaseVisit) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 4.2)
bxp <- ggboxplot(mydata, x = "DiseaseVisit", y = "shannon", color = "Disease")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
dev.off()


png(filename=paste0(filename,"shannon.DiseaseVisit.wilcox.png"), width=800, height = 800)
stat.test <- mydata %>%
  #group_by(Country) %>%
  wilcox_test(shannon ~ DiseaseVisit) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 4.2)
bxp <- ggboxplot(mydata, x = "DiseaseVisit", y = "shannon", color = "Disease")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
dev.off()

wilcox.shannon <- pairwise.wilcox.test(mydata$shannon, 
                                       mydata$DiseaseVisit, 
                                       p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("shannon -- DiseaseVisit wilcoxon adjustp")
tab.shannon

write.table(tab.shannon, file = "MulSanne.shannon.wilcox.p.txt", sep = "\t", row.names = FALSE)





################################################################
KP_Disease=kruskal.test(observed ~ Disease, data = mydata)$p.value
KP_Visit=kruskal.test(observed ~ Visit, data = mydata)$p.value
KP_DiseaseVisit=kruskal.test(observed ~ DiseaseVisit, data = mydata)$p.value
KP_Sid=kruskal.test(observed ~SID, data = mydata)$p.value

print(paste0("KP_Disease=", KP_Disease, 
             "; KP_Visit=", KP_Visit,
             "; KP_Sid=", KP_Sid, 
             "; KP_DiseaseVisit=", KP_DiseaseVisit
     ))

#all data : "KP_Disease=2.01e-06; KP_Visit=0.32; KP_Sid=4.58e-05; KP_DiseaseVisit=0.000124"
################################################################
KP_Disease=kruskal.test(shannon ~ Disease, data = mydata)$p.value
KP_Visit=kruskal.test(shannon ~ Visit, data = mydata)$p.value
KP_DiseaseVisit=kruskal.test(shannon ~ DiseaseVisit, data = mydata)$p.value
KP_Sid=kruskal.test(shannon ~SID, data = mydata)$p.value

print(paste0("KP_Disease=", KP_Disease, 
             "; KP_Visit=", KP_Visit,
             "; KP_Sid=", KP_Sid, 
             "; KP_DiseaseVisit=", KP_DiseaseVisit
))
#all data: [1] "KP_Disease=0.000195 KP_Visit=0.14; KP_Sid=0.00107; KP_DiseaseVisit=0.00328"


x<-test_d+1 ###to avoid 0 
beta_dist <- vegdist(x, index = "bray")
mds <- metaMDS(beta_dist)
mds_data <- cbind(as.data.frame(mds$points), mydata)
mds_data$SampleID <- rownames(mds_data)

#install.packages('devtools')
library(devtools)
#install_github('fawda123/ggord')
library(ggord)
library(ggplot2)
# https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html


png(filename=paste0(filename,".beta.DiseaseVisit.png"),  width=600, height=380)
p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = DiseaseVisit)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=DiseaseVisit, fill=DiseaseVisit), alpha=0.1, geom="polygon")
print(p1)
dev.off() 

png(filename=paste0(filename,".beta.Disease.Visit.png"),  width=600, height=380)
p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Disease)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=Disease, fill=Disease), alpha=0.1, geom="polygon")+ facet_grid(. ~ Visit)
print(p1)
dev.off() 

png(filename=paste0(filename,".beta.Visit.Disease.png"),  width=600, height=380)
p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Visit)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=Visit, fill=Visit), alpha=0.1, geom="polygon")+ facet_grid(. ~ Disease)
print(p1)
dev.off() 

png(filename=paste0(filename,".beta.Disease.SID.png"),  width=580, height=3000)
p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Disease)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=Disease, fill=Disease), alpha=0.1, geom="polygon")+ facet_grid(SID ~ Disease)
print(p1)
dev.off() 

#png(filename=paste0(filename,".beta.Visit.SID.png"),  width=580, height=3000)
#p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Visit)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=Visit, fill=Visit), alpha=0.1, geom="polygon")+ facet_grid(SID ~ Visit )
#print(p1)
#dev.off() 

png(filename=paste0(filename,".beta.SID.png"),  width=580, height=3000)
p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = SID)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=SID, fill=SID), alpha=0.1, geom="polygon")+ facet_grid(SID ~ . )
print(p1)
dev.off() 

############################################################


#print("all samples")
# calculate Bray-Curtis distance among samples
#test_d2 <- test_d[BodySite !="Control",]
bc.dist <- vegdist(test_d, method = "bray", na.rm = TRUE)  # This removes missing values
#bc.dist <- vegdist(test_d2, method = "bray")
# cluster communities using average-linkage algorithm
bc.clust <- hclust(bc.dist, method = "average")

# plot cluster diagram
png(filename=paste0(filename,".bray-Curtis-cluster.png"), width=2000, height=600)
plot(bc.clust, ylab = "Bray-Curtis dissimilarity")
dev.off()
# Taxonomic (Bray-Curtis) dissimilarity explained
#print("adonis result bc.dist ~ DiseaseVisit")
adonis2(bc.dist ~ DiseaseVisit) #p=0.001
adonis2(bc.dist ~ Disease) #p=0.003
adonis2(bc.dist ~ Visit) #p=0.001


######################################
jaccard.dist <- vegdist(x, method = "jaccard")
mds2 <- metaMDS(jaccard.dist)
mds_data2 <- as.data.frame(cbind(mds2$points, mydata))
mds_data2$SampleID <- rownames(mds_data2)

png(filename=paste0(filename,".beta.DiseaseVisit.jaccard.png"),  width=380, height=380)
p1<-ggplot(mds_data2, aes(x = MDS1, y = MDS2, color = DiseaseVisit)) +geom_point(size=6)+theme_bw()+
  stat_chull(aes(color=DiseaseVisit, fill=DiseaseVisit), alpha=0.1, geom="polygon") 
print(p1)
dev.off()
png(filename=paste0(filename,".beta.Disease.Visit.jaccard.png"),  width=800, height=380)
p1<-ggplot(mds_data2, aes(x = MDS1, y = MDS2, color = Disease)) +geom_point(size=6)+theme_bw()+
  stat_chull(aes(color=Disease, fill=Disease), alpha=0.1, geom="polygon") + facet_grid(. ~ Visit)
print(p1)
dev.off()
png(filename=paste0(filename,".beta.Visit.Disease.jaccard.png"),  width=800, height=380)
p1<-ggplot(mds_data2, aes(x = MDS1, y = MDS2, color = Visit)) +geom_point(size=6)+theme_bw()+
  stat_chull(aes(color=Visit, fill=Visit), alpha=0.1, geom="polygon") + facet_grid(. ~ Disease)
print(p1)
dev.off()

adonis2(jaccard.dist ~ DiseaseVisit) #p=0.001
adonis2(jaccard.dist ~ Disease) #p=0.001
adonis2(jaccard.dist ~ Visit) #p=0.001

# cluster communities using average-linkage algorithm
#bc.clust <- hclust(bc.dist, method = "average")
# plot cluster diagram
#png(filename=paste0(filename,".bray-Curtis-cluster.png"), width=2000, height=800)
#plot(bc.clust, ylab = "Bray-Curtis dissimilarity")
#dev.off()


#test.adonis <- adonis2(bc.dist ~ DiseaseVisit)
#test.adonis <- as.data.frame(test.adonis$aov.tab)
#test.adonis
#p=0.001

###################################################################
#PAIRWISE PERMANOVA
cbn <- combn(x=unique(DiseaseVisit), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- x [c(DiseaseVisit %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(DiseaseVisit %in% cbn[,i]), ]
  
  permanova_pairwise <- adonis2(vegdist(ps.subs, method = "bray") ~ metadata_sub$DiseaseVisit)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
# Output the results to a file
write.table(p.table, "beta_DiseaseVisit.txt", sep = "\t", quote = FALSE)











###################################################################
#PAIRWISE PERMANOVA
cbn <- combn(x=unique(DiseaseVisit), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- x [c(metadata_sub$DiseaseVisit %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(metadata_sub$DiseaseVisit %in% cbn[,i]), ]
  
  permanova_pairwise <- adonis2(vegdist(ps.subs, method = "bray") ~ metadata_sub$DiseaseVisit)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
# Output the results to a file
write.table(p.table, "beta_DiseaseVisit.txt", sep = "\t", quote = FALSE)


###########################Jaccard test
#print("adonis result jaccard.dist ~ RashPH")
#adonis2(jaccard.dist ~ DiseaseVisit)
#test.adonis <- adonis2(jaccard.dist ~ DiseaseVisit)



