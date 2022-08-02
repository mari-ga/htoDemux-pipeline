install.packages("Matrix", repos="http://R-Forge.R-project.org")
install.packages("ggpubr")
## Load Seurat library
library(Seurat)
## Load other libraries needed fr Seurat library(dplyr)
library(Matrix)
library(ggpubr)
library(gplots)
library(argparser, quietly=TRUE)


data_dir <- '/home/icb/mariana.gonzales/storage/MS_nuclei_hashing/747495_GX12/outs/filtered_feature_bc_matrix'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
seurat_object = CreateSeuratObject(counts = expression_matrix)


