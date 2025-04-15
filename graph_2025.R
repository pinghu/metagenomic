####################
##Ping Hu 4-15-2025
######################


rm(list=ls())
library(ggplot2)
library(stringr)
#filename="metaphlan.relab10.7"
filename="gladiator_zpt_short"

library("ggpubr")
library(ggtext)
library(xfun)
library(rstatix)
library(dplyr)
library(tidyverse)

library(openxlsx) 

A<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A);
C=A[1:d[1], 2:d[2]]

Cname=colnames(A)[2:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
Treat=rep("NA", Clen)
Sid=rep("NA", Clen)
Visit=rep("NA", Clen)
TreatVisit=rep("NA", Clen)
ASFS=rep("NA", Clen)
age=rep("NA", Clen)
id2=rep("NA", Clen)
splitname<-strsplit(Cname, "[.]")
for(mm in  1:Clen ){
  Treat[mm]=splitname[[mm]][2]
  if(splitname[[mm]][2] =="A"){
    Treat[mm]="A:ZPT"
  }else if(splitname[[mm]][2] =="B"){
    Treat[mm]="B:PO"
  }else if(splitname[[mm]][2] =="C"){
    Treat[mm]="C:Ctl"
  }
  Visit[mm]=splitname[[mm]][3]
  Sid[mm]=splitname[[mm]][4]
  ASFS[mm]=as.numeric(splitname[[mm]][5])
  age[mm]=splitname[[mm]][6]
  id2[mm]=splitname[[mm]][7]
  TreatVisit[mm]=paste0(Treat[mm], Visit[mm])
}
ASFS=as.numeric(ASFS[ASFS !="NA"])

library("plotrix")
truefc<-function(VVV){
  #print(VVV)
  if (is.finite(VVV )){
    XXX=VVV
    if(VVV==0){
      XXX=NA
    }else if(VVV<1){
      XXX=-1/VVV
    }
    return(XXX)
  }else{
    return("NA")
  }
}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  
  summary_func <- function(x, col){
    vals <- x[[col]]
    vals_clean <- vals[!is.na(vals)]
    c(
      mean = mean(vals_clean),
      sd = sd(vals_clean),
      se = sd(vals_clean) / sqrt(length(vals_clean)),
      n = length(vals_clean)
    )
  }
  
  data_sum <- ddply(data, groupnames, .fun = summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  
  return(data_sum)
}


trt_groups=unique(TreatVisit)
# # Generate all possible pairwise comparisons
trt_pairs <- combn(trt_groups, 2, simplify = FALSE)

big_results_list <- list()

for (i in 1:d[1]){
  
  genename=A[i,1]
  gene<-as.numeric(C[i,])
  splitG<-strsplit(as.character(genename), "[.]")
  LLL=length(splitG[[1]])
  genus=splitG[[1]][LLL]
  cleaned_genus <- gsub("^\\w__", "", genus)  # Removes the prefix
  cleaned_genus <- gsub("_", " ", cleaned_genus)  # Replaces underscores with spaces
  # ccc=as.numeric(min(gene[gene>0&!is.na(gene)]))/10
  # RA=gene+ccc
  # logRA=log(RA)
  mydata=data.frame(gene, Treat, Visit, Sid, ASFS, id2,TreatVisit)

  df2 <- mydata %>%
    group_by(Treat, Visit) %>%
    summarize(
      mean_RA = mean(gene, na.rm = TRUE),
      se_RA = sd(gene, na.rm = TRUE) / sqrt(n()),
      count=n()
    )
  
  
  
  df2$TreatVisit=paste0(df2$Treat, df2$Visit)
  
  KP_TreatVisit=kruskal.test(gene ~ TreatVisit, data = mydata)$p.value
  percentZeroNA=(sum(C[i,]==0)+sum(is.na(C[i,])))/(d[2]-1)
  meanAll <-mean(gene)
  results_list <- list()
  results_list[["TaxonName"]] <-cleaned_genus
  results_list[["TaxonFullName"]] <-genename
  results_list[["PercentZero"]] <-percentZeroNA
  results_list[["MeanAll"]] <-meanAll
  
  for ( test in trt_groups) { 
    results_list[[paste0("mean.", test)]] <- df2$mean_RA[df2$TreatVisit == test]
    results_list[[paste0("se.", test)]] <- df2$se_RA[df2$TreatVisit == test]
    results_list[[paste0("sample_number.", test)]] <- df2$count[df2$TreatVisit == test]
    
  }
  results_list[["Kruskal_P"]] <-KP_TreatVisit
  
  
  for (pair in trt_pairs) {
    
    trt1 <- pair[1]
    trt2 <- pair[2]
    
    # Extract the gene values for each treatment group
    data1 <- mydata[mydata$TreatVisit == trt1, "gene"]
    data2 <- mydata[mydata$TreatVisit == trt2, "gene"]
    
    if (length(data1) > 1 && length(data2) > 1) {
      
      m1 <- mean(data1, na.rm = TRUE)
      m2 <- mean(data2, na.rm = TRUE)
      fc <- m1 / m2
      
      TFC<- truefc(fc)
      
      # Check if there's variation before running Wilcoxon
      if (length(unique(data1)) > 1 && length(unique(data2)) > 1) {
        wilcox <- wilcox.test(data1, data2, paired = FALSE)$p.value
      } else {
        wilcox <- NA
      }
      
      # Store results with named keys
      label <- paste0(trt1, "_vs_", trt2)
      results_list[[paste0("wilcox.", label)]] <- wilcox
      results_list[[paste0("truefold.", label)]] <- TFC
    }
  }
  
  
  big_results_list[[genename]] <- results_list
  
  
  
  df2$TreatOriginal=df2$Treat
  df2$Treat[df2$Treat == "A:ZPT"] <- "ZPT (N=50)"
  df2$Treat[df2$Treat == "C:Ctl"] <- "Placebo (N=50)"
  

  
 
  # Map colors to each group manually
  treat_colors <- c("ZPT (N=50)" = "#1f77b4",   # Blue
                    "Placebo (N=50)" = "#ff7f0e") # Orange
  
  # Add colored strip labels to your dataframe
  df2$Treat_label <- ifelse(df2$Treat == "ZPT (N=50)",
                            "<span style='color:#1f77b4'>ZPT (N=50)</span>",
                            "<span style='color:#ff7f0e'>Placebo (N=50)</span>")
  
  png(filename = paste0("kraken.", i, ".", genus, ".png"), width = 600, height = 480, res = 180)
  
  p <- ggplot(df2, aes(x = Visit, y = mean_RA, group = Treat, color = Treat)) + 
    geom_line() +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_RA - se_RA, ymax = mean_RA + se_RA), width = 0.2,
                  position = position_dodge(0.05)) +
    facet_wrap(~Treat_label) +
    labs(x = "Visit", y = paste0(genus, "(%)")) +
    scale_color_manual(values = treat_colors) +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_markdown(face = "bold"),
      axis.title.y = element_text(face = "italic")  # Italicize y-axis label
    )
  
  print(p)
  dev.off()
  
  
  
  # png(filename=paste0("kraken.",i,".",genus, ".pair.png"),  width=600, height=480, res=180)
  # p1=ggplot(data=mydata,aes(x=Visit,y=gene,color=Treat, group=Sid)) +ylab(paste(genus, "Abundance"))+geom_point(size=1)+geom_line(color="grey")+ theme_classic()+theme(legend.position="none")+facet_wrap(~Treat)
  # 
  # print(p1)
  # 
  # dev.off()
  # ########## I need to add the statistical significant into this
  # png(filename = paste0(genus, ".png"), width = 1000, height = 1000, res = 300)
  # bxp <- ggboxplot(mydata, x = "TreatVisit", y = "gene", color = "Treat") +
  #   ggtitle(paste0(cleaned_genus,"\nKP=", sprintf("%.02f", KP_TreatVisit), " 0%=", sprintf("%.02f", 100*percentZeroNA))) + ####This is for pathway
  #   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
  #   ylab(paste(cleaned_genus, "(%)")) + 
  #   theme(
  #     legend.position = "none",
  #     axis.text.x = element_text(angle = 60, hjust = 1),
  #     plot.title = element_text(face = "italic"), # Italicizes the entire title
  #     axis.title.y = element_text(face = "italic")  # Italicize y-axis label
  #   )
  # print(bxp)
  # dev.off()
  
  
}
# Convert big results list to a single data frame
big_results_table <- bind_rows(big_results_list)
pData <- big_results_table %>% select("Kruskal_P", starts_with("wilcox"))
pData[is.na(pData)] <- 1
d=dim(pData)
df<-data.frame(matrix(ncol=d[2], nrow=d[1]))
colnames(df)=paste0("fdr.",colnames(pData))

for (i in 1:d[2]) {
  p <- as.numeric(pData[[i]])       # Use [[i]] or as.numeric() to avoid list issues
  df[, i] <- p.adjust(p, method = "fdr")
}
final_result_table <- cbind(big_results_table, df)
write.table(final_result_table, file = paste0(filename, ".stat"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)     
#library(openxlsx)  # Load the package
output_file <- paste0(filename, ".stat.xlsx")
write.xlsx(final_result_table, file = output_file, rowNames = FALSE, colNames = TRUE)

print(paste("Results successfully written to:", output_file))