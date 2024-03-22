##############################################################
###http://t-redactyl.io/blog/2016/02/creating-plots-in-r-using-ggplot2-part-6-weighted-scatterplots.html
###################################################################

rm(list = ls())
library(ggplot2)
library(dplyr)
#file="Porphyromonas_gingivalis_T00145.pstat.MP"

args <- commandArgs(trailingOnly = TRUE)
print(args)
file <- args[1]
#file="map00550.koenrich"
A<-read.table(file,header=TRUE, sep="\t",stringsAsFactors=F)
d=dim(A)
colnames(A)[2]="Treat"
colnames(A)[1]="Strain"
A=A[A$ID!="ID",]
A=A[A$Strain != "Porphyromonas_gingivalis_W83",]
A=A[A$Strain != "Porphyromonas_gingivalis__ATCC33277",]
A$minuslogp= -log10(as.numeric(A$pvalue))

png(filename = paste0(file, ".bubble.png"), width = 1800, height = 800, res = 120)

# Generate the plot
ggplot(A, aes(y = Strain, x = Treat, size = minuslogp, color = Direction)) +
  geom_point(aes(shape = Direction), stroke = 2) +  # Map shape to Direction
  scale_shape_manual(values = c("All" = 16, "Up" = 24, "Down" = 25)) +  # Define shapes for each Direction
  scale_size_area(max_size = 15) +  # Adjust the size scale
  ggtitle(paste(A$ID[1], A$Description[1])) +  # Set the title based on the first entries in A$ID and A$Description
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "All" = "orange")) +  # Map colors
  theme(legend.position = "bottom", legend.direction = "horizontal",  # Adjust legend position
        axis.title.x = element_text(size = 14, face = "bold"),  # Make x-axis title bold and larger
        axis.title.y = element_text(size = 14, face = "bold")) +  # Make y-axis title bold and larger
  theme_bw()  # Use the base theme for a cleaner look

# Close the device used for file output
dev.off()




