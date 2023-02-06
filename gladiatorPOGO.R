# ---
#   title: "Gladiator_gene_graph"
# author: "Ping Hu"
# date: "2022-1-30"
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
#  rm(args)

filename ="GladiatorPO.go.stat.xls"
B<-read.table(filename, sep="\t", header=TRUE)
dB <- dim(B);
dB
A=filter(B, B$fdr.B_PCor_Spearman <=0.1)
d <- dim(A);
d
x=filter(B, B$fdr.C_PCor_Spearman <=0.1)
d <- dim(x);
d
A1=filter(A, A$fdr.Pwilcox_Boct <=0.1)
d <- dim(A1);
d
A2=filter(A1, A1$fdr.Pwilcox_BvC_log <=0.1)
d <- dim(A1);
d

A1=filter(A2, abs(A2$truefc_B) >=2)
d <- dim(A1);
d
A2=filter(A1, abs(A1$MaxAve) >=0.00009)
d <- dim(A2);
d

write.csv(A2, file=paste0(filename, ".selected022623.csv"))

#A=filter(B, B$fdr.B_PCor_Spearman <=0.05
            # & B$fdr.PCor_Spearman <=0.05
            # & B$fdr.B_PCor_Pearson <=0.05
#             & (B$fdr.Pwilcox_Boct <=0.05)
             #& (abs(B$CorASFS_Spearman)>=0.3)
#             & (abs(B$B_CorASFS_Spearman)>=0.4)
#             & (abs(B$truefc_B)>=2)
#             & (B$fdr.Pwilcox_BvC <=0.05)
#             & (B$MaxAve>=0.00009)
#             & grepl("[BP]", B$TaxonName)
#)

A=filter(A2, grepl("[BP]", A2$GOCategory))
d <- dim(A);
d
#library(simplifyEnrichment)
#go_id=A$GOID
#length(go_id)
#mat = GO_similarity(go_id, ont='BP')
#png(filename = paste0(filename,".go.cluster.png"), width = 3380, height = 1860, res=300)
#df = simplifyGO(mat, word_cloud_grob_param = list(max_width = 80))
#dev.off()

library(rrvgo)
simMatrix <- calculateSimMatrix(A$GOID,
                                #orgdb="org.Hs.eg.db",
                                orgdb = "org.Sc.sgd.db",
                                #orgdb="org.EcK12.eg.db",
                                ont="BP",
                                method="Rel")
#scores <- setNames(-log10(A$Pwilcox_BvC_log), A$GOID)
#scores <- setNames(-log10(A$fdr.Pwilcox_Boct), A$GOID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                #scores,
                                threshold=0.7,
                                #orgdb="org.Hs.eg.db"
                                orgdb="org.Sc.sgd.db"
                                #orgdb="org.EcK12.eg.db"
                                )
dim(reducedTerms)

part1=reducedTerms[reducedTerms$go== reducedTerms$parent,]
dim(part1)
part1$term

part1$parent=""
part1$parentTerm=""
part2=reducedTerms[reducedTerms$go!= reducedTerms$parent,]
dim(part1)
newdf=rbind(part1, part2)

heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)
png(filename = paste0(filename,".revigoScatter.png"), width = 2660, height = 1860, res=300)
scatterPlot(simMatrix, reducedTerms)
dev.off()
png(filename = paste0(filename,".revigoTreePlot.png"), width = 2660, height = 1860, res=300)
treemapPlot(reducedTerms)
dev.off()
wordcloudPlot(reducedTerms, min.freq=1, colors="black")

write.csv(reducedTerms, file=paste0(filename, ".reducedTerms022623.csv"))
#library(webr)
#PieDonut(reducedTerms,aes(pies=parent,donuts=go))
library(plotly)
fig <- plot_ly()


#if (!require("processx")) install.packages("processx")
fig <- fig %>% add_trace(
  type='sunburst',
  ids=newdf$term,
  labels=newdf$term,
  parents=newdf$parentTerm,
  domain=list(column=1),
  maxdepth=2,
  insidetextorientation='radial'
)
#orca(fig,paste0(filename,".sunburst.svg") )
#plotly::export(p = fig, #the graph to export
#               paste0(filename,".sunburst.png")) #the name and type of file (can be .png, .jpeg, etc.)fig
library(htmlwidgets)
htmlwidgets::saveWidget(
  widget = fig, #the plotly object
  file = "sunburst.html", #the path & file name
  selfcontained = TRUE #creates a single html file
)
A=filter(A, A$GOID %in% reducedTerms$go)
dim(A)
joined_df <- merge(A, reducedTerms, by.x = "GOID", 
                   by.y = "go", all.x = TRUE, all.y = FALSE)
A=joined_df
p1_table= bind_cols( #transmute(A, taxon = A$GOName) ,
                     transmute(A, taxon = paste0(A$parent, "--", A$GOName)) ,
                     transmute(A, fold_or_cor = truefc_B), 
                     transmute(A, diff_cat = fdr.Pwilcox_Boct), 
                     transmute(A, diff_p = Pwilcox_Boct)
                     )%>%
  mutate(group = "1:PO W3vsBL")

p2_table= bind_cols(  #transmute(A, taxon = A$GOName) ,
                     transmute(A, taxon = paste0(A$parent, "--", A$GOName)) ,
                     transmute(A, fold_or_cor = truefc_C),
                     transmute(A, diff_cat = fdr.Pwilcox_Cctl), 
                     transmute(A, diff_p = Pwilcox_Cctl)
                     )  %>%
  mutate(group = "3:Ctl W3vsBL")
p3_table=bind_cols(  #transmute(A, taxon = A$GOName) ,
                    transmute(A, taxon = paste0(A$parent, "--", A$GOName)) ,
                    transmute(A, fold_or_cor = B_CorASFS_Spearman),
                    transmute(A, diff_cat = fdr.B_PCor_Spearman), 
                    transmute(A, diff_p = B_PCor_Spearman)
                    ) %>%
  mutate(group = "2:PO ASFS Correlation")

p4_table=bind_cols(  #transmute(A, taxon = A$GOName) ,
                    transmute(A, taxon = paste0(A$parent, "--", A$GOName)) ,
                    transmute(A, fold_or_cor = C_CorASFS_Spearman),
                    transmute(A, diff_cat = fdr.C_PCor_Spearman),
                    transmute(A, diff_p = C_PCor_Spearman) 
                    ) %>%
  mutate(group = "4:Ctl ASFS Correlation")

PAll=bind_rows(p1_table,p3_table,p2_table,p4_table)
#PAll=bind_rows(p1_table,p3_table)
PAll$diff_cat[is.na(PAll$diff_cat)]=1
PAll$diff_cat[is.na(PAll$diff_p)]=1
PAll$diff_Cor=recode((PAll$diff_p <=0.05)*1, 
                     `1` = "Pvalue <=0.05",
                     `0` = "Pvalue >0.05, Not significant")

PAll$diff_Cor=recode((PAll$diff_cat <=0.1)*1, 
                     `1` = "fdr Pvalue <0.1, Significant",
                     `0` = PAll$diff_Cor)
PAll$fold_or_cor[PAll$fold_or_cor < -10]=-10 
PAll$fold_or_cor[PAll$fold_or_cor > 10]=10 
PAll=PAll[order(PAll$fold_or_cor),]
png(filename = paste0(filename,".asfs_fdr_reduced.png"), width = 4880, height = 2860, res=300)
p<-    ggplot(PAll,aes(y = taxon, x = fold_or_cor, color = diff_Cor)) +geom_point(size=3) +
  geom_vline(xintercept = 0, color = "yellow", size=1)+
  scale_color_brewer(palette = "Dark2", name = NULL) +
  facet_grid(cols  = vars(group), scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+ 
  theme(legend.position = "bottom")
print(p)
dev.off()

######################Not Very Useful, term too general#########################

###############################################################
simMatrix <- calculateSimMatrix(A$GOID,
                                #orgdb="org.Hs.eg.db",
                                #orgdb = "org.Sc.sgd.db",
                                orgdb="org.EcK12.eg.db",
                                ont="BP",
                                method="Rel")
#scores <- setNames(-log10(A$Pwilcox_BvC_log), A$GOID)
#scores <- setNames(-log10(A$fdr.Pwilcox_Boct), A$GOID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                #scores,
                                threshold=0.7,
                                #orgdb="org.Hs.eg.db"
                                #orgdb="org.Sc.sgd.db"
                                orgdb="org.EcK12.eg.db"
)
dim(reducedTerms)

part1=reducedTerms[reducedTerms$go== reducedTerms$parent,]
dim(part1)
part1$term

part1$parent=""
part1$parentTerm=""
part2=reducedTerms[reducedTerms$go!= reducedTerms$parent,]
dim(part1)
newdf=rbind(part1, part2)

heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)
png(filename = paste0(filename,".revigoScatter.bac.png"), width = 2660, height = 1860, res=300)
scatterPlot(simMatrix, reducedTerms)
dev.off()
png(filename = paste0(filename,".revigoTreePlot.bac.png"), width = 2660, height = 1860, res=300)
treemapPlot(reducedTerms)
dev.off()
wordcloudPlot(reducedTerms, min.freq=1, colors="black")

write.csv(reducedTerms, file=paste0(filename, ".reducedTerms022623.bac.csv"))
#library(webr)
#PieDonut(reducedTerms,aes(pies=parent,donuts=go))
library(plotly)
fig <- plot_ly()


#if (!require("processx")) install.packages("processx")
fig <- fig %>% add_trace(
  type='sunburst',
  ids=newdf$term,
  labels=newdf$term,
  parents=newdf$parentTerm,
  domain=list(column=1),
  maxdepth=2,
  insidetextorientation='radial'
)
#orca(fig,paste0(filename,".sunburst.svg") )
#plotly::export(p = fig, #the graph to export
#               paste0(filename,".sunburst.png")) #the name and type of file (can be .png, .jpeg, etc.)fig
library(htmlwidgets)
htmlwidgets::saveWidget(
  widget = fig, #the plotly object
  file = "sunburst.bac.html", #the path & file name
  selfcontained = TRUE #creates a single html file
)
