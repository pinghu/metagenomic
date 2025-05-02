##################################################################
# Ping Hu
# 4-17-2025 Direct generate the powerpoint slides to save time
#################################################################
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(tidyverse)
library("plotrix")
library(vegan)
library(rlang)
#install.packages("officer")
#install.packages("rvg")         # Optional: for editable vector plots
library(officer)
library(rvg)
#install.packages('devtools')
library(devtools)
#install_github('fawda123/ggord')
library(ggord)
library(flextable)
library(dendextend)
#library(officer)
#filename1= "GSS3121_meta_data_4_18_2025.txt"
filename1="GSS3121_allPassed_meta.txt"
args <- commandArgs(trailingOnly = TRUE)
#filename<- args[1]
filename="metaphlan.estimated_count.filtered.7"
A2<-read.table(filename1, sep="\t", header=TRUE)
d <- dim(A2);
order_levels <- c("Shave", "Trim")
A2$TestTime <- factor(A2$Trimeatment, levels = order_levels, ordered = TRUE)


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

test <- read.table(filename, header = TRUE, row.names = 1, sep="\t")
test[is.na(test)]<-0
d=dim(test)
d
X=apply(test, 1, sum)
test_filter2=test[(X/d[2] >=10),]

test_d=data.frame(t(test))


min(apply(test_d, 1, sum)) ###162338

test_a <- decostand(test_d, method = "total")

Cname=rownames(test_d)
Cname_clean <- sub("\\.metaphlan4$", "", Cname)
matching_dataset <- A2 %>%
  filter(GenomicID %in% Cname_clean) %>%
  arrange(match(GenomicID, Cname_clean))
# Step 1: Remove .metaphlan4 from rownames of test_d
rownames(test_d) <- sub("\\.metaphlan4$", "", rownames(test_d))
# Step 2: Create a named vector for lookup from matching_dataset
# Ensure both columns exist
id_map <- setNames(matching_dataset$NewID, matching_dataset$GenomicID)
# Step 3: Replace rownames using the mapping
# Only replace if the name exists in the mapping
updated_rownames <- id_map[rownames(test_d)]
# Step 4: Keep existing name if no match was found (fallback to original name)
rownames(test_d) <- ifelse(!is.na(updated_rownames), updated_rownames, rownames(test_d))
#https://stackoverflow.com/questions/30057765/histogram-ggplot-show-count-label-for-each-bin-for-each-category
shannon<-diversity(test_d)
simp<-diversity(test_d, "simpson")
invsimp<-diversity(test_d, "inv")

observed<-apply(test_d>0,1,sum)
N <- apply(test_d,1,sum)
###Richness Index#####
Menhinick_index<-observed/sqrt(N)
Margalef_index <-(observed-1)/log(N)

mydata0=data.frame(shannon,simp,invsimp, observed, Menhinick_index, Margalef_index ,Cname_clean)
mydata <- inner_join(mydata0, matching_dataset, by = c("Cname_clean" = "GenomicID"))


write.table(t(mydata), file =paste0(filename, ".alphadiversity.RShort"), col.names=TRUE, row.names=TRUE, sep = "\t")
x<-test_d
beta_dist <- vegdist(x, index = "bray")
mds <- metaMDS(beta_dist)
mds_data <- cbind(as.data.frame(mds$points), mydata)
mds_data$SampleID <- rownames(mds_data)
# calculate Bray-Curtis distance among samples
bc.dist <- vegdist(test_d, method = "bray")
# cluster communities using average-linkage algorithm
bc.clust <- hclust(bc.dist, method = "average")

dend <- as.dendrogram(bc.clust)
# Step 3: Match sample order with TestTime groups
sample_order <- labels(dend)
dendrogram_png <- paste0(filename, ".bray-Curtis-cluster.png")
png(filename = dendrogram_png, width = 3000, height = 1600, res = 300)
# Set larger bottom margin for label space: mar = c(bottom, left, top, right)
par(mar = c(12, 4, 4, 2))  # increase bottom space
# Plot dendrogram with larger label space and readable size
plot(dend,
     main = "Bray-Curtis Cluster",
     ylab = "Dissimilarity",
     cex = 1.0,                  # control label text size
     las = 2                     # rotate labels (optional)
)
dev.off()

doc <- read_pptx()
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "Cluster Dendrogram", 
          location = ph_location_type(type = "title")) %>%
  ph_with(external_img(dendrogram_png, width = 9, height = 5), 
          location = ph_location_type(type = "body"))

# Add slide with both plots and caption
# Loop over all columns in mydata
for (colname in colnames(mydata)) {
  
  if (colname == "shannon") next
  if (colname == "observed") next
  # Check how many unique values
  # Convert the variable to a factor explicitly in the plot data
  
  
  n_unique <- length(unique(mydata[[colname]]))
  if (n_unique < 2) next
  # Proceed only if unique values <= 10
  if (n_unique <= 10) {
    message("Plotting for grouping variable: ", colname)
    mydata[[colname]] <- as.factor(mydata[[colname]])
    # Use sym() to capture the column as a symbol
    xvar <- sym(colname)
    
    body_site_factor <- as.factor(mydata[[colname]])[match(sample_order, mydata$NewID)]
    
    # Assign group-based label colors
    group_palette <- rainbow(length(levels(body_site_factor)))
    labels_colors(dend) <- group_palette[body_site_factor]
    
    # Optional: Set label text (if you want to rename labels)
    # labels(dend) <- mydata$NewID[match(labels(dend), mydata$NewID)]
    
    # Step 4: Plot nicely with more space
    dendrogram_png <- paste0(filename, ".bray-Curtis-cluster.",colname,".png")
    png(filename = dendrogram_png, width = 3000, height = 1600, res = 300)
    
    # Set larger margin for long sample names
    par(mar = c(12, 4, 4, 2))
    
    # Plot
    plot(dend,
         main = "Bray-Curtis Cluster by Body Site",
         ylab = "Dissimilarity",
         cex = 1.0,      # label size
         las = 2         # rotate label text
    )
    
    # Add legend
    legend("topright", legend = levels(body_site_factor), fill = group_palette, border = NA, bty = "n")
    
    dev.off()
    
    doc <- doc %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with(value = paste0("Cluster Dendrogram colored by ", colname),
              location = ph_location_type(type = "title")) %>%
      ph_with(external_img(dendrogram_png, width = 9, height = 5), 
              location = ph_location_type(type = "body"))
    
    
    
    
    
    
    # Open PNG device
    png_filename <- paste0(filename, ".shannon.", colname, ".png")
    png(filename = png_filename, width = 1600, height = 1600, res = 300)


    stat.test <- tryCatch({
       mydata %>%
        t_test(as.formula(paste("shannon ~", colname))) %>%
        mutate(
          y.position = 1.05 * max(mydata$shannon, na.rm = TRUE),
          label = ifelse(p <= 0.05, paste0("p = ", signif(p, 2)), "")
          )
    }, error = function(e) {
      message("Skipping stat.test for ", colname, ": ", e$message)
      return(NULL)
    })
    # Base plot
    p <- ggplot(mydata, aes(x = !!xvar, y = shannon, color = !!xvar)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.5) +
      geom_jitter(width = 0.2, height = 0, size = 2.5, alpha = 0.7) +
      labs(
        title = paste("Shannon Diversity by", colname),
        y = "Shannon Diversity", x = colname
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12, face = "bold")
      )
    if (!is.null(stat.test) && nrow(stat.test) > 0) {
      p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, step.increase=0.1)
    }
    
    print(p)
    dev.off()
    
    png_filename <- paste0(filename, ".observed.", colname, ".png")
    png(filename = png_filename, width = 1600, height = 1600, res = 300)
    # Run t-test and filter for significant p-values
   
    stat.test <- tryCatch({
      mydata %>%
        t_test(as.formula(paste("observed ~", colname))) %>%
        mutate(
          y.position = 1.01 * max(mydata$observed, na.rm = TRUE),
          label = ifelse(p <= 0.05, paste0("p = ", signif(p, 2)), "")
          )
    }, error = function(e) {
      message("Skipping stat.test for ", colname, ": ", e$message)
      return(NULL)
    })
    
    # Base plot
    p2 <- ggplot(mydata, aes(x = !!xvar, y = observed, color = !!xvar)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.5) +
      geom_jitter(width = 0.2, height = 0, size = 2.5, alpha = 0.7) +
      labs(
        title = paste("Observed Species by", colname),
        y = "Observed Species", x = colname
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12, face = "bold")
      )
    if (!is.null(stat.test) && nrow(stat.test) > 0) {
      p2 <- p2 + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, step.increase=0.1)
    }
    print(p2)
    dev.off()

    
    KP_shannon=kruskal.test(as.formula(paste("shannon ~", colname)), data = mydata)$p.value 
    KP_observed=kruskal.test(as.formula(paste("observed ~", colname)) , data = mydata)$p.value 
    caption_text <- paste0("Shannon Diversity kruskal wallis test pvalue: ",colname,"=", sprintf("%.2f", KP_shannon), 
                           "; Observed Species kruskal wallis test p value:",colname,"=", sprintf("%.2f", KP_observed))
    print(caption_text)
    doc <- doc %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with(value = "Alpha Diversity Summary", location = ph_location_type(type = "title")) %>%
      # Add left plot (bxp1)
      ph_with(dml(ggobj = p), location = ph_location(left = 0.5, top = 1.2, width = 4.5, height = 4.5)) %>%
      # Add right plot (bxp2)
      ph_with(dml(ggobj = p2), location = ph_location(left = 5.2, top = 1.2, width = 4.5, height = 4.5)) %>%
      # Add caption text with smaller font
      ph_with(
        fpar(ftext(caption_text, prop = fp_text(font.size = 12))),
        location = ph_location(left = 0.5, top = 6.2, width = 9, height = 1)
      )
    
    #####Now beta diversity####
    # Filter out rows with NA in the grouping variable
    non_na_idx <- which(!is.na(mydata[[colname]]))
    # Subset the distance matrix and metadata
    bc.dist.sub <- as.dist(as.matrix(bc.dist)[non_na_idx, non_na_idx])
    mydata.sub <- mydata[non_na_idx, , drop = FALSE]
    mds_data.sub <- mds_data[non_na_idx, , drop = FALSE]
    # Make sure there is more than one group with at least two samples
    group_counts <- table(mydata.sub[[colname]])
    if (length(group_counts) < 2 || any(group_counts < 2)) {
      message("Skipping ", colname, ": not enough valid group sizes for adonis")
    } else {
      png(filename=paste0(filename,".beta.",colname,".png"),  width=1800, height=800, res=300)
      p3<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = !!xvar )) +geom_point(size=3, alpha = 0.7)+
        theme_bw()+stat_chull(aes(color= !!xvar, fill= !!xvar), alpha=0.1, geom="polygon")
      print(p3)
      dev.off() 
      
      # Build the formula as a string
      group_formula <- as.formula(paste("bc.dist.sub ~", colname))
      # Run adonis
      adonis_result <- adonis2(group_formula, data = mydata.sub)
      # Extract stats
      pval <- adonis_result$`Pr(>F)`[1]
      sample_size <- nrow(mydata.sub)
      # Format summary
      adonis_text <- paste0(colname, " Bray-Curtis PERMANOVA: n = ", sample_size,
                            " R² = ", round(adonis_result$R2[1], 3),
                            ", p = ", format.pval(pval, digits = 3, eps = 0.001))
     # print(adonis_result)  # Or use in caption
      print(adonis_text)  # Or use in caption
      #PAIRWISE PERMANOVA
      cbn <- combn(x=unique(mds_data.sub[[colname]]), m = 2)
      p <- c()
      for(i in 1:ncol(cbn)){
        ps.subs <- x [c(mds_data[[colname]] %in% cbn[,i]), ]
        metadata_sub <- mds_data[c(mds_data[[colname]] %in% cbn[,i]), ]
        permanova_pairwise <- adonis2(vegdist(ps.subs, method = "bray") ~ metadata_sub[[colname]])
        p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
      }
      if (!is.null(p)) {
      p.adj <- p.adjust(p, method = "BH")
      p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
      p.table
      # Output the results to a file
      write.table(p.table, paste0(filename, ".beta.", colname, ".txt"), sep = "\t", quote = FALSE)
      
      # Create a formatted flextable
      pairwise_ft <- flextable(p.table) %>%
        autofit() %>%
        set_table_properties(width = 1, layout = "autofit") %>%
        theme_booktabs() %>%
        align(align = "center", part = "all") %>%
        fontsize(size = 9, part = "all")
      }else{
        pairwise_ft=""
      }
      doc <- doc %>%
        add_slide(layout = "Title and Content", master = "Office Theme") %>%
        
        # Custom title with font size 36
        ph_with(
          value = fpar(ftext("Beta Diversity with Bray-Curtis Distance & Pairwise Permanova Test",
                             fp_text(font.size = 36))),
          location = ph_location(left = 0.5, top = 0.3, width = 9, height = 1)
        ) %>%
        
        # Add the ggplot object
        ph_with(dml(ggobj = p3), 
                location = ph_location(left = 0.5, top = 2.0, width = 5.5, height = 4.2)) %>%
        
        # Add the caption text under the plot
        ph_with(fpar(ftext(adonis_text, prop = fp_text(font.size = 16))), 
                location = ph_location(left = 0.5, top = 6.3, width = 5.5, height = 0.5)) %>%
        
        # Add the pairwise permanova table next to the plot
        ph_with(pairwise_ft, 
                location = ph_location(left = 6.2, top = 2.0, width = 3.5, height = 5))
      
    }
    
  }
}

print(doc, target = "diversity_summary.pptx")
