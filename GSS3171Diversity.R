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

#install.packages("officer")
#install.packages("rvg")         # Optional: for editable vector plots

library(officer)
library(rvg)

filename1= "GSS3171.meta_merge.clean.txt"
args <- commandArgs(trailingOnly = TRUE)
#filename<- args[1]
filename="metaphlan.estimated_count.filtered.7.clean"
A2<-read.table(filename1, sep="\t", header=TRUE)
d <- dim(A2);
order_levels <- c("Day1AM",  "Day1PM",  "Day2",    "Day8")
A2$TestTime <- factor(A2$Time, levels = order_levels, ordered = TRUE)


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

#saveRDS(mydata,paste0(filename, ".diversitydata"))
#saveRDS(test_d,paste0(filename,".test_d"))
#https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/
png(filename=paste0(filename, ".shannon.Time.png"), width=1600, height=1600, res=300)
stat.test1 <- mydata %>%
  t_test(shannon ~ TestTime) %>%
  mutate(y.position = 1.05*max(mydata$shannon))
bxp1 <- ggboxplot(mydata, x = "TestTime", y = "shannon", color = "TestTime")
bxp1<- bxp1 + stat_pvalue_manual(stat.test1,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  ggtitle("shannon Diversity : t-test")+
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) #+ scale_x_discrete(limits = order_levels)
print(bxp1)
dev.off()

########################observed#################################
png(filename=paste0(filename, ".observed.Treat.png"), width=1600, height=2000, res=300)
stat.test <- mydata %>%
  t_test(observed ~ TestTime) %>%
  mutate(y.position = 1.05*max(mydata$observed))
bxp2 <- ggboxplot(mydata, x = "TestTime", y = "observed", color = "TestTime")
bxp2<- bxp2 + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  ggtitle("observed Species: t-test")+
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) #+ scale_x_discrete(limits = order_levels)
print(bxp2)
dev.off()

######################################################33333
KP_TestTime_shannon=kruskal.test(shannon ~ TestTime, data = mydata)$p.value 
KP_TestTime_observed=kruskal.test(observed ~ TestTime, data = mydata)$p.value 
caption_text <- paste0("Shannon Diversity kruskal wallis test pvalue: TestTime=", sprintf("%.2f", KP_TestTime_shannon), 
                       "; Observed Species kruskal wallis test p value: TestTime=", sprintf("%.2f", KP_TestTime_observed))

print(caption_text)
#[1] "Shannon KP_Treat=0.00135089426057143; KP_Treat_observed=0.010213571184125"
# Create PowerPoint
doc <- read_pptx()

# Add slide with both plots and caption
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "Alpha Diversity Summary", location = ph_location_type(type = "title")) %>%
  
  # Add left plot (bxp1)
  ph_with(dml(ggobj = bxp1), location = ph_location(left = 0.5, top = 1.2, width = 4.5, height = 4.5)) %>%
  
  # Add right plot (bxp2)
  ph_with(dml(ggobj = bxp2), location = ph_location(left = 5.2, top = 1.2, width = 4.5, height = 4.5)) %>%
  
  # Add caption text with smaller font
  ph_with(
    fpar(ftext(caption_text, prop = fp_text(font.size = 12))),
    location = ph_location(left = 0.5, top = 6.2, width = 9, height = 1)
  )

########################################
x<-test_d
beta_dist <- vegdist(x, index = "bray")
mds <- metaMDS(beta_dist)
mds_data <- cbind(as.data.frame(mds$points), mydata)
mds_data$SampleID <- rownames(mds_data)

#install.packages('devtools')
library(devtools)
#install_github('fawda123/ggord')
library(ggord)
library(ggplot2)

 png(filename=paste0(filename,".beta.TestTime.png"),  width=1800, height=800, res=300)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = TestTime)) +geom_point(size=3, alpha = 0.7)+
   theme_bw()+stat_chull(aes(color=TestTime, fill=TestTime), alpha=0.1, geom="polygon")
 print(p1)
 dev.off() 

 
############################################################
print("all samples")
# calculate Bray-Curtis distance among samples
bc.dist <- vegdist(test_d, method = "bray")
# cluster communities using average-linkage algorithm
bc.clust <- hclust(bc.dist, method = "average")
# plot cluster diagram

# Run the adonis2 test
adonis_result <- adonis2(bc.dist ~ mydata$TestTime)
# Extract p-value
pval <- adonis_result$`Pr(>F)`[1]
# Get sample size (non-NA rows in bc.dist and TestTime)
sample_size <- sum(!is.na(mydata$TestTime))
# Build summary string
adonis_text <- paste0("Bray-Curtis PERMANOVA: n = ", sample_size,
                      ", R² = ", round(adonis_result$R2[1], 3),
                      ", p = ", format.pval(pval, digits = 3, eps = 0.001))

library(dendextend)
# Step 1: Hierarchical clustering
bc.clust <- hclust(bc.dist, method = "average")  # or method = "ward.D2"
# Step 2: Convert to dendrogram object
dend <- as.dendrogram(bc.clust)
# Step 3: Get grouping (TestTime) in correct order
# Match TestTime to the order of leaves
sample_order <- labels(dend)
group_colors <- as.factor(mydata$TestTime)[match(sample_order, rownames(mydata))]
# Assign colors to groups (customize as needed)
group_palette <- rainbow(length(levels(group_colors)))  # or use scale_brewer, etc.
labels_colors(dend) <- group_palette[group_colors]
# Optional: Color branches
dend <- color_branches(dend, k = length(levels(group_colors)), groupLabels = levels(group_colors))
# # Step 4: Plot the dendrogram
# png(filename = paste0(filename, ".bray-Curtis-cluster.png"), width = 3000, height = 1600, res = 300)
# plot(dend, main = "Bray-Curtis Cluster", ylab = "Dissimilarity")
# legend("topright", legend = levels(group_colors), fill = group_palette, border = NA, bty = "n")
# dev.off()

# Save the plot as a PNG
dendrogram_png <- paste0(filename, ".colored_dendrogram.png")
png(filename = dendrogram_png, width = 3000, height = 1600, res = 300)
plot(dend, main = "Colored Bray-Curtis Cluster", ylab = "Dissimilarity")
legend("topright", legend = levels(group_colors), fill = group_palette, bty = "n")
dev.off()

doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "Cluster Dendrogram Colored by TestTime", 
          location = ph_location_type(type = "title")) %>%
  ph_with(external_img(dendrogram_png, width = 9, height = 5), 
          location = ph_location_type(type = "body"))


# ######################################
# jaccard.dist <- vegdist(x, method = "jaccard")
# mds2 <- metaMDS(jaccard.dist)
# mds_data2 <- as.data.frame(cbind(mds2$points, mydata))
# mds_data2$SampleID <- rownames(mds_data2)
# 
# png(filename=paste0(filename,".beta.PH_Buttocks_Category.jaccard.png"),  width=1200, height=800, res=300)
# p1<-ggplot(mds_data2, aes(x = MDS1, y = MDS2, color = PH_Buttocks_Category)) +
#   geom_point(size=3, alpha=0.7)+theme_bw()+
#   stat_chull(aes(color=PH_Buttocks_Category, fill=PH_Buttocks_Category), alpha=0.1, geom="polygon") 
# print(p1)
# dev.off()
# 
# # test.adonis <- adonis2(bc.dist ~ RashPH)
# # test.adonis <- as.data.frame(test.adonis$aov.tab)
# # test.adonis
###################################################################
#PAIRWISE PERMANOVA
cbn <- combn(x=unique(mds_data$TestTime), m = 2)
p <- c()
for(i in 1:ncol(cbn)){
  ps.subs <- x [c(mds_data$TestTime %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(mds_data$TestTime %in% cbn[,i]), ]
  permanova_pairwise <- adonis(vegdist(ps.subs, method = "bray") ~ metadata_sub$TestTime)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
# Output the results to a file
write.table(p.table, paste0(filename, ".beta_Treat.txt"), sep = "\t", quote = FALSE)

library(flextable)
#library(officer)

# Example: assume you already have pairwise_table
# Example columns: Group1 | Group2 | R2 | p.value | p.adjust

# Create a formatted flextable
pairwise_ft <- flextable(p.table) %>%
  autofit() %>%
  set_table_properties(width = 1, layout = "autofit") %>%
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  fontsize(size = 9, part = "all")


# Create the flextable if not already done
pairwise_ft <- flextable(p.table) %>%
  autofit() %>%
  set_table_properties(width = 1, layout = "autofit") %>%
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  fontsize(size = 9, part = "all")


doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  
  # Custom title with font size 36
  ph_with(
    value = fpar(ftext("Beta Diversity with Bray-Curtis Distance & Pairwise Permanova Test",
                       fp_text(font.size = 36))),
    location = ph_location(left = 0.5, top = 0.3, width = 9, height = 1)
  ) %>%
  
  # Add the ggplot object
  ph_with(dml(ggobj = p1), 
          location = ph_location(left = 0.5, top = 2.0, width = 5.5, height = 4.2)) %>%
  
  # Add the caption text under the plot
  ph_with(fpar(ftext(adonis_text, prop = fp_text(font.size = 16))), 
          location = ph_location(left = 0.5, top = 6.3, width = 5.5, height = 0.5)) %>%
  
  # Add the pairwise permanova table next to the plot
  ph_with(pairwise_ft, 
          location = ph_location(left = 6.2, top = 2.0, width = 3.5, height = 5))

print(doc, target = "diversity_summary.pptx")
