rm(list=ls())
# Install and load the required packages
#install.packages(c("readxl", "pheatmap"))
library(readxl)
library(pheatmap)
####Here are different color choices#######https://github.com/WishartLab/heatmapper/blob/master/expression/server.R
library("RColorBrewer")
lowCol = "#0016DB"
midCol = "#FFFFFF"
highCol = "#FFFF00"
num_colors = 75
BYcol=colorRampPalette(c(lowCol, midCol, highCol))(num_colors)
lowCol = "#FF0000"
midCol = "#000000"
highCol = "#33FF00"
RGcol=colorRampPalette(c(lowCol, midCol, highCol))(num_colors)
lowCol = "#000000"
midCol = "#bdbdbd"
highCol = "#FFFFFF"
Graycol=colorRampPalette(c(lowCol, midCol, highCol))(num_colors)
lowCol = "#C9438C"
midCol = "#f7f7f7"
highCol = "#7BC134"
Piygcol=colorRampPalette(c(lowCol, midCol, highCol))(num_colors)
#rainbow_col=rainbow(num_colors)
#topo_col=topo.colors(num_colors)
#####################################################################
# Load the data from the Excel file
args <- commandArgs(trailingOnly = TRUE)
print(args)
file_path <- args[1]
outname <- args[2]
rm(args)
#file_path <- "test.xlsx"  # Replace with the actual path to your Excel, column with Name as title is used for row name
data <- read_excel(file_path)

# Remove "PWY-" prefix from genenames for better display
#data$genename <- gsub("PWY-", "", data$genename)

# Select columns for heatmap
heatmap_data <- data.matrix(data[,-1]);
row.names(heatmap_data)<-data$Name
# Identify numeric columns
#numeric_cols <- sapply(heatmap_data, is.numeric)

# Convert data to a matrix
heatmap_matrix <- as.matrix(heatmap_data)

# Row normalization
heatmap_matrix_normalized <- t(scale(t(heatmap_matrix)))

# Perform clustering
#clustering_rows <- hclust(dist(heatmap_matrix_normalized))
#clustering_cols <- hclust(dist(t(heatmap_matrix_normalized)))

# Create a heatmap using pheatmap
png(filename=paste0(outname, "pheat.png"), height=3000, width=2600, res=200)
pheatmap(heatmap_matrix_normalized,  color = BYcol, main = "Heatmap of Pathways")
dev.off()

png(filename=paste0(outname, "pheat.cut.png"), height=3000, width=2600, res=200)
pheatmap(heatmap_matrix_normalized,  color = BYcol,cutree_rows = 6,
         cutree_cols = 2, main = "Heatmap of Pathways")
dev.off()

#the following normalization is same as above
#cal_z_score <- function(x){
#  (x - mean(x)) / sd(x)
#}

#data_subset_norm <- t(apply(heatmap_matrix, 1, cal_z_score))

#library("amap")
#library("pvclust")
#library("gplots")


#Not as nice as pheatmap
#png(filename=paste0("heat.png"), height=3000, width=2600, res=200)
#par(mar = c(5, 5, 4, 20))
#par(mai = c(1.5, 1.5, 1.5, 1.5))
#heatmap.2(heatmap_matrix,trace='none', density='none', scale='row', labRow=data$Name,col=BYcol,margins=c(15,12), cexRow = 0.5)
#dev.off()

#my_hclust_gene <- hclust(dist(heatmap_matrix), method = "complete")

# install if necessary
#install.packages("dendextend")

# load package
#library(dendextend)
#png(filename=paste0("2.dendex.png"), height=3000, width=2600, res=200)
#as.dendrogram(my_hclust_gene) %>%
#  plot(horiz = TRUE)
#dev.off()
#my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 6)