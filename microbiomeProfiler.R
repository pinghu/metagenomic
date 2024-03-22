rm(list=ls())
library(MicrobiomeProfiler)
#filename="test.KO"
args <- commandArgs(trailingOnly = TRUE)
print(args)
filename <- args[1]
A<-read.table(filename, sep="\t", header=TRUE)
FC <-as.numeric(A[,2])
P <-as.numeric(A[,3])

sig <- A[P <=0.05,1]
sigUp <- A[P<=0.05&FC>0,1]
sigDown <- A[FC<0&P<=0.05,1]

combined_results <- data.frame()
x <- enrichKO(sig)
df_x =data.frame(x)
if(nrow(df_x)>0){
  df_x$Direction="All"
  combined_results <-rbind(combined_results, df_x)
}
y <- enrichKO(sigUp)
df_y=data.frame(y)
if(nrow(df_y)>0){
  df_y$Direction="Up"
  combined_results <-rbind(combined_results, df_y)
}
z <- enrichKO(sigDown)
df_z=data.frame(z)
if(nrow(df_z)>0){
  df_z$Direction="Down"
  combined_results <-rbind(combined_results, df_z)
}

write.table(combined_results, file = paste0(filename, ".MFko.xls"), sep = "\t", row.names = FALSE, col.names = TRUE, eol = "\n", na = "NA")

