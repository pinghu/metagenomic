# ---
#   title: "Gladiator_gene_graph"
# author: "Ping Hu"
# date: "2022-11-16"
# output: html_document
#####Note: taxon clean not clean , need to automate no matter which level it is if 6 need to remove species
  ###prepare data -- matrix for correlation analysis
  #############
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
#print(args)
filename <- args[1]
#  outname <- args[2]
#  rm(args)
#filename ="vox1.metaphlan3.count10.7.all.stat.xls"
#filename ="vox1.metaphlan3.count10.7.bacteria.stat.xls"
#filename ="vox1.kraken.count10.7.fungal.stat.xls"
#filename ="vox1.pathcoverage.community.stat.xls"
filename ="vox1.ec.stat.ann.xls"
B<-read.table(filename, sep="\t", header=TRUE)
d <- dim(B);
A=filter(B,B$Pwilcox_DN <=0.05 |  B$Pval_ASFS_COR_spearman <0.05 )
#A=filter(A,A$aveAll >=0.01  )
d <- dim(A);
write.csv(A, file=paste0(filename, ".sigP.csv"))


p_table = bind_cols(
  transmute(A, taxon = A[,1]) ,
  transmute(A, fold_or_cor = truefc_DN) ,
  transmute(A, diff_cat = Pwilcox_DN)) %>%
  mutate(group = "Dandruff_vs_NoDandruff") 



p4_table=bind_cols( transmute(A, taxon = A[,1]) ,transmute(A, fold_or_cor = ASFS_Cor_Spearman),
                    transmute(A, diff_cat = Pval_ASFS_COR_spearman)) %>%
  mutate(group = "ASFS Correlation")

PAll=bind_rows(p_table,p4_table)
PAll$diff_cat[is.na(PAll$diff_cat)]=1

PAll$diff_Cor=recode((PAll$diff_cat <0.05)*1, 
                     `1` = "Significant",
                     `0` = "Nonsignificant")
PAll$fold_or_cor[PAll$fold_or_cor < -10]=-10 
PAll$fold_or_cor[PAll$fold_or_cor > 10]=10 




png(filename = paste0(filename,".fold_cor_sig.png"), width = 2880, height = 4880, res=300)
p<-    ggplot(PAll,aes(y = taxon, x = fold_or_cor, color = diff_Cor)) +geom_point() +
  geom_vline(xintercept = 0, color = "yellow", size=1)+
  scale_color_brewer(palette = "Dark2", name = NULL) +
  facet_grid(cols  = vars(group), scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(p)
dev.off()
