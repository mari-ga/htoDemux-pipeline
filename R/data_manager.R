install.packages("Matrix", repos="http://R-Forge.R-project.org")
install.packages("ggpubr")
## Load Seurat library
library(Seurat)
## Load other libraries needed fr Seurat library(dplyr)
library(Matrix)
library(ggpubr)
library(gplots)


data_dir_umi <- '/home/icb/mariana.gonzales/storage/MS_nuclei_hashing/747495_GX12/outs/filtered_feature_bc_matrix'

list.files(data_dir_umi) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
umi <- Read10X(data.dir = data_dir_umi)
seurat_object = CreateSeuratObject(counts = umi)
saveRDS(seurat_object, file="/home/icb/mariana.gonzales/storage/tranformed_data/GX12_UMI.rds")


data_dir_counts <- '/home/icb/mariana.gonzales/storage/MS_nuclei_hashing/combined_GX34/outs/filtered_feature_bc_matrix'
list.files(data_dir_counts) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
counts <- Read10X(data.dir = data_dir_counts)
seurat_object_2 = CreateSeuratObject(counts = counts)
saveRDS(seurat_object_2, file="/home/icb/mariana.gonzales/storage/tranformed_data/GX12_counts.rds")




