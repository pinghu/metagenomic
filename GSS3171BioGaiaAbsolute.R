##########################
#Ping Hu April 16, 2025
##############################
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
filename2 = "GSS3171.meta_merge.clean.txt"
A2<-read.table(filename2, sep="\t", header=TRUE)
d2 <- dim(A2);
filename <- args[1]
#filename="metaphlan.relab10.6.clean"
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(tidyverse)
library(ggtext)
library(xfun)

library(openxlsx) 
library(circlize)
library(ComplexHeatmap)
library(lattice)   
library(reshape2)
library(gridExtra)

library(ggplot2)
library(stringr)


my.t.test.p.value <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}

my.wilcox.p.value <- function(...) {
    obj<-try(wilcox.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}
my.kruskal.p.value <- function(...) {
  obj<-try(kruskal.test(...), silent=TRUE)
  if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}

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

A0<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A0);
B=A0[1:d[1], 2:d[2]]
ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
C=B+ZZ

Cname=colnames(A0)[2:d[2]]
Cname_clean <- sub("\\.metaphlan4$", "", Cname)
order_levels <- c("Day1AM",  "Day1PM",  "Day2",    "Day8")
A2$TestTime <- factor(A2$Time, levels = order_levels, ordered = TRUE)
matching_dataset <- A2 %>%
  filter(GenomicID %in% Cname_clean) %>%
  arrange(match(GenomicID, Cname_clean))


trt_groups=unique(A2$TestTime)
# # Generate all possible pairwise comparisons
trt_pairs <- combn(trt_groups, 2, simplify = FALSE)
big_results_list <- list()

MicroRNAYield <- as.numeric(A2$BacterialDNA_ng)
names(MicroRNAYield) <- A2$GenomicID
# Match values using OrgID in matching_dataset
matching_dataset$MicroRNAYield <- MicroRNAYield[match(matching_dataset$GenomicID, names(MicroRNAYield))]
for (i in 1:d[1]){
  genename=A0[i,1]
  gene_relative <-as.numeric(C[i,])
  meanAllRelative <-mean(gene_relative)
  splitG<-strsplit(as.character(genename), "[|]")
  LLL=length(splitG[[1]])
  genus=splitG[[1]][LLL]
  cleaned_genus <- gsub("^\\w__", "", genus)  # Removes the prefix
  cleaned_genus <- gsub("_", " ", cleaned_genus)  # Replaces underscores with spaces
  percentZeroNA=(sum(B[i,]==0)+sum(is.na(B[i,])))/(d[2]-1)
  mydata0=data.frame(gene_relative, Cname_clean)
  mydata <- inner_join(mydata0, matching_dataset, by = c("Cname_clean" = "GenomicID"))
  mydata$gene_absolute=mydata$gene_relative * mydata$MicroRNAYield /100
  meanAllAbsolute <-mean(mydata$gene_absolute)
  KP_Treat_relative=kruskal.test(gene_relative ~ TestTime, data = mydata)$p.value
  KP_Treat_absolute=kruskal.test(gene_absolute ~ TestTime, data = mydata)$p.value
  
  summary_relative <- mydata %>%
    group_by(TestTime) %>%
    summarize(
      mean_RA = mean(gene_relative, na.rm = TRUE),
      se_RA = sd(gene_relative, na.rm = TRUE) / sqrt(n()),
      count=n()
    )
  
  summary_absolute <- mydata %>%
    group_by(TestTime) %>%
    summarize(
      mean_RA = mean(gene_absolute, na.rm = TRUE),
      se_RA = sd(gene_absolute, na.rm = TRUE) / sqrt(n()),
      count=n()
    )
  
  # Initialize results list for pairwise tests
  results_list <- list()
  results_list[["TaxonName"]] <-cleaned_genus
  results_list[["TaxonFullName"]] <-genename
  results_list[["PercentZero"]] <-percentZeroNA
  results_list[["MeanAllRelative"]] <-meanAllRelative
  results_list[["MeanAllAbsolute_ng"]] <-meanAllAbsolute
  
  for ( test in trt_groups) { 
    results_list[[paste0("mean_relative.", test)]] <- summary_relative$mean_RA[summary_relative$TestTime == test]
    results_list[[paste0("se_relative.", test)]] <- summary_relative$se_RA[summary_relative$TestTime == test]
    results_list[[paste0("sampleN_relative.", test)]] <- summary_relative$count[summary_relative$TestTime == test]
    results_list[[paste0("mean_absolute.", test)]] <- summary_absolute$mean_RA[summary_absolute$TestTime == test]
    results_list[[paste0("se_absolute.", test)]] <- summary_absolute$se_RA[summary_absolute$TestTime == test]
    results_list[[paste0("sampleN_absolute.", test)]] <- summary_absolute$count[summary_absolute$TestTime == test]
  }
  results_list[["Kruskal_P.relative"]] <-KP_Treat_relative
  results_list[["Kruskal_P.absolute"]] <-KP_Treat_absolute
  
  for (pair in trt_pairs) {
    trt1 <- pair[1]
    trt2 <- pair[2]
    # Extract the gene values for each treatment group
    data1 <- mydata[mydata$TestTime == trt1, "gene_relative"]
    data2 <- mydata[mydata$TestTime == trt2, "gene_relative"]
    data1A <- mydata[mydata$TestTime == trt1, "gene_absolute"]
    data2A <- mydata[mydata$TestTime == trt2, "gene_absolute"]
    if (length(data1) > 1 && length(data2) > 1) {
      m1 <- mean(data1, na.rm = TRUE)
      m2 <- mean(data2, na.rm = TRUE)
      m1A <- mean(data1A, na.rm = TRUE)
      m2A <- mean(data2A, na.rm = TRUE)     
      fc_relative <- m1 / m2
      TFC_relative <- truefc(fc_relative)
      fc_absolute <- m1A / m2A
      TFC_absolute <- truefc(fc_absolute)      
      if (length(unique(data1)) > 1 && length(unique(data2)) > 1) {
        wilcox_relative <- wilcox.test(data1, data2, paired = FALSE)$p.value
      } else {
        wilcox_relative <- 1
      }
      if (length(unique(data1A)) > 1 && length(unique(data2A)) > 1) {
        wilcox_absolute <- wilcox.test(data1A, data2A, paired = FALSE)$p.value
      } else {
        wilcox_absolute <- 1
      }
      # Store results with named keys
      label <- paste0(trt1, "_vs_", trt2)
      results_list[[paste0("wilcox_relative.", label)]] <- wilcox_relative
      results_list[[paste0("truefold_relative.", label)]] <- TFC_relative
      results_list[[paste0("wilcox_absolute.", label)]] <- wilcox_absolute
      results_list[[paste0("truefold_absolute.", label)]] <- TFC_absolute      
    }
  }
  
  # Store the transformed results for this gene
  big_results_list[[genename]] <- results_list
  
  p <- ggplot(summary_absolute, aes(x = TestTime, y = mean_RA, group = TestTime, color = TestTime)) + 
    geom_line(group=1) +
    geom_point(size = 3) +
    ggtitle(paste0(cleaned_genus," ng \nKP=", sprintf("%.02f", KP_Treat_absolute), " 0%=", sprintf("%.02f", 100*percentZeroNA))) +
    geom_errorbar(aes(ymin = mean_RA - se_RA, ymax = mean_RA + se_RA), width = 0.2,
                  position = position_dodge(0.05)) +
    labs(x = "Visit", y = paste0(cleaned_genus, "(ng):absolute")) +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_markdown(face = "bold"),
      axis.title.y = element_text(face = "italic")  # Italicize y-axis label
    )
  
  p2 <- ggplot(summary_relative, aes(x = TestTime, y = mean_RA, group = TestTime, color = TestTime)) + 
    geom_line(group=1) +
    geom_point(size = 3) +
    ggtitle(paste0(cleaned_genus," %\nKP=", sprintf("%.02f", KP_Treat_relative), " 0%=", sprintf("%.02f", 100*percentZeroNA))) + 
    geom_errorbar(aes(ymin = mean_RA - se_RA, ymax = mean_RA + se_RA), width = 0.2,
                  position = position_dodge(0.05)) +
    labs(x = "Visit", y = paste0(cleaned_genus, "(%):relative")) +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_markdown(face = "bold"),
      axis.title.y = element_text(face = "italic")  # Italicize y-axis label
    ) 
  
  png(filename = paste0(genus, ".combine.png"), width = 2000, height = 1000, res = 300)
  print(grid.arrange(p, p2, ncol = 2))
  dev.off()

      
}


# Convert big results list to a single data frame
big_results_table <- bind_rows(big_results_list)


pData <- big_results_table %>% select(starts_with("Kruskal"), starts_with("wilcox"))
pData[is.na(pData)] <- 1
d=dim(pData)
df<-data.frame(matrix(ncol=d[2], nrow=d[1]))
colnames(df)=paste0("fdr.",colnames(pData))

for (i in 1:d[2]) {
  p <- as.numeric(pData[[i]])       # Use [[i]] or as.numeric() to avoid list issues
  df[, i] <- p.adjust(p, method = "fdr")
}
final_result_table <- cbind(big_results_table, df)
write.table(final_result_table, file = paste0(filename, ".absolute_relavent_stat"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)     

output_file <- paste0(filename, ".absolute_relavent.stat.xlsx")
write.xlsx(final_result_table, file = output_file, rowNames = FALSE, colNames = TRUE)

print(paste("Results successfully written to:", output_file))


# Filter the big_results_table for rows where KP_Treat_AR <= 0.05
filtered_results <- big_results_table[(big_results_table$Kruskal_P.relative <= 0.05 | big_results_table$Kruskal_P.absolute <=0.05),]
output_file2 <- paste0(filename, ".", dim(filtered_results)[1], ".sig.xlsx")
write.xlsx(filtered_results, file = output_file2, rowNames = FALSE, colNames = TRUE)

meanRelative <-final_result_table  %>% select( starts_with("mean_relative"))
rownames(meanRelative) =final_result_table$TaxonName
meanRelative_zscore_mat <- t(apply(meanRelative, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))

# Remove rows with any NA, NaN, or Inf
meanRelative_zscore_mat_clean <- meanRelative_zscore_mat[complete.cases(meanRelative_zscore_mat) & 
                                                 apply(meanRelative_zscore_mat, 1, function(x) all(is.finite(x))), ]

# Now plot
htcorRelative <- Heatmap(meanRelative_zscore_mat_clean, name = "Zscore Relative Abundance", rect_gp = gpar(col = "white", lwd = 2),
                    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
                    row_names_side = "left",
                    column_names_gp = gpar(fontsize = 8),
                    column_title = "Relative Percentage",
                    show_row_dend = FALSE,
                    show_column_dend = FALSE,
                    cluster_columns = FALSE
)

meanAbsolute <-final_result_table  %>% select( starts_with("mean_absolute"))
rownames(meanAbsolute) =final_result_table$TaxonName
meanAbsolute_zscore_mat <- t(apply(meanAbsolute, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))

# Remove rows with any NA, NaN, or Inf have to use relavent filter to make show row number is consistent
meanAbsolute_zscore_mat_clean <- meanAbsolute_zscore_mat[complete.cases(meanRelative_zscore_mat) & 
                                                           apply(meanRelative_zscore_mat, 1, function(x) all(is.finite(x))), ]

meanAbsolute_clean <-meanAbsolute[complete.cases(meanRelative_zscore_mat) & 
                                    apply(meanRelative_zscore_mat, 1, function(x) all(is.finite(x))), ]

meanRelative_clean <-meanRelative[complete.cases(meanRelative_zscore_mat) & 
                                    apply(meanRelative_zscore_mat, 1, function(x) all(is.finite(x))), ]

# Now plot
htcorAbsolute <- Heatmap(meanAbsolute_zscore_mat_clean, name = "Zscore Absolute ng", rect_gp = gpar(col = "white", lwd = 2),
                         row_names_gp = gpar(fontsize = 8, fontface = "italic"),
                         row_names_side = "left",
                         column_names_gp = gpar(fontsize = 8),
                         column_title = "Absolute ng",
                         show_row_dend = FALSE,
                         show_column_dend = FALSE,
                         cluster_columns = FALSE
)

htAbsolute <- Heatmap(meanAbsolute_clean, name = "Absolute ng", rect_gp = gpar(col = "white", lwd = 2),
                         row_names_gp = gpar(fontsize = 8, fontface = "italic"),
                         row_names_side = "left",
                         column_names_gp = gpar(fontsize = 8),
                         column_title = "Absolute ng",
                         show_row_dend = FALSE,
                         show_column_dend = FALSE,
                         cluster_columns = FALSE
)

htRelative <- Heatmap(meanRelative_clean, name = "Relative %", rect_gp = gpar(col = "white", lwd = 2),
                      row_names_gp = gpar(fontsize = 8, fontface = "italic"),
                      row_names_side = "left",
                      column_names_gp = gpar(fontsize = 8),
                      column_title = "Relative %",
                      show_row_dend = FALSE,
                      show_column_dend = FALSE,
                      cluster_columns = FALSE
)
ht_list=htcorRelative+htcorAbsolute +htRelative+ htAbsolute
jpeg(paste0(filename, ".absolute.jpg"), width = 2500, height = 1200, res=150 )
draw(ht_list, heatmap_legend_side = "right", 
     column_title = paste( filename, "mean_ZScore"))
dev.off()
########################## comparison summary figures#########################
meanData <-filtered_results  %>% select( starts_with("mean_"))
rownames(meanData) =filtered_results$TaxonName

fcData <- filtered_results %>% select(starts_with("truefold"))
rownames(fcData) =filtered_results$TaxonName
# Convert all columns to numeric
#fcData[] <- lapply(fcData, function(col) as.numeric(as.character(col)))
newfcData <- as.matrix(fcData)
mode(newfcData) <- "numeric" 

pData <- filtered_results %>% select( "TaxonFullName", starts_with("Kruskal"), starts_with("wilcox"))
pData[is.na(pData)] <- 1


newfcData[newfcData > 5] <- 5
newfcData[newfcData < -5] <- -5

pData <- filtered_results %>% select("TaxonFullName", starts_with("wilcox_AR"))
adjData <-filtered_results %>% select("TaxonFullName", starts_with("adj"))
create_pk_mat <- function(pk, adjpk) {
  nrow_pk <- nrow(pk)
  ncol_pk <- ncol(pk)
  pk_mat <- matrix("", nrow = nrow_pk, ncol = ncol_pk)
  # Set the row and column names
  rownames(pk_mat) <- rownames(pk)
  colnames(pk_mat) <- colnames(pk)
  pk_mat[!is.na(pk) & pk <= 0.1] <- "*"
  pk_mat[!is.na(pk) & pk <= 0.05] <- "**"
  pk_mat[!is.na(adjpk) & adjpk <= 0.1] <- "***"
  pk_mat[!is.na(adjpk) & adjpk <= 0.05] <- "****"
  return(pk_mat)
}
pmat <-create_pk_mat(pData, adjData)

ht4 <- Heatmap(fcData,name="Fold Change", 
               column_names_gp = gpar(fontsize = 9), 
               cluster_columns = FALSE, 
               row_names_side = "left", 
               rect_gp = gpar(col = "white", lwd = 2),
               row_names_gp = gpar(fontsize = 8, fontface = "italic"), 
               cluster_rows = FALSE, 
               cell_fun = function(j, i, x, y, w, h, col) {
                 grid.text(pmat[i, j], x, y)
               },
               col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")))

png(paste0(filename, ".foldheat.png"), width = 1200, height = 1200, res = 150)  # Adjust size and resolution as needed
#print(htcor)
draw(ht4)
# Close the PNG device
dev.off()

