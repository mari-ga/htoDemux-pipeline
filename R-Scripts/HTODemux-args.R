#!/usr/bin/env Rscript
install.packages("Seurat",repos=c("http://cran.rstudio.com"))
install.packages('multtest')
#Receive arguments from command line
options(echo=TRUE)
#Get arguments as a vector
myargs = commandArgs(trailingOnly=TRUE)

#Library
library(Seurat)

#Import files
file_umis = myargs[0]
print("This is the file:")
file_umis

pbmc.umis <-readRDS(file_umis)
print(pbmc.umis)

