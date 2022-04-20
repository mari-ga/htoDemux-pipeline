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
file_htos = myargs[2]



#---------------- Section 1 - Input files -----------------
pbmc.umis <-readRDS(file_umis)
print(pbmc.umis)

pbmc.htos <- readRDS(file_htos)
#pbmc.htos

#Identify which UMI corresponds to which hashtag.
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Confirm that the HTO have the correct names
rownames(pbmc.htos)
