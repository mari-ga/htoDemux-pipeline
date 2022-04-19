#!/usr/bin/env Rscript
install.packages("Seurat",repos=("http://cran.rstudio.com"))
install.packages("spatstat.sparse",repos=("http://cran.rstudio.com"))
#Receive arguments from command line
options(echo=TRUE)
#Get arguments as a vector
myargs = commandArgs(trailingOnly=TRUE)

#Library
library(Seurat)

#Import files
file_umis = myargs[1]
typeof(file_umis)

pbmc.umis <-readRDS(file_umis)
print(pbmc.umis)

