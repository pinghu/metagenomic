##############################################################
#This R code is intended to renormalize the bacterial number #
#into relative abundance which will added up to 1 columnwise #
#Ping HU 06/03/2016
#remmeber need toshift the first row as the name
##############################################################

#filename="pathcoverage.community.tsv"
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outname <- args[2]

counts <- read.table(filename, header=TRUE, sep='\t', comment.char='', quote='', row.names=1)
counts <- as.data.frame(t(counts))
counts <- counts[order(rownames(counts)), ]

# to look at part of your dataframe, use head() or counts[1:3,1:3] Divide by row sums to get relative abundances.sweep(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...),x is array, margin dimension ( 2=column,1=row),

abundances  <- sweep(counts, 1, rowSums(counts), FUN='/')

# Inspect,you can use head or specify dimension.

abundances <- as.data.frame(t(abundances))
write.table(abundances, file=paste0(filename, ".norm"), col.names=TRUE, sep="\t", row.names=TRUE)
