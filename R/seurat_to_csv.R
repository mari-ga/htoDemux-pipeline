library(Seurat)
library(data.table)



counts <- Read10X(data.dir = "/home/icb/mariana.gonzales/storage/MS_nuclei_hashing/combined_GX34/outs/filtered_feature_bc_matrix")


data_to_write_out <- as.data.frame(as.matrix(counts))
fwrite(x = data_to_write_out, row.names = TRUE, file = "/home/icb/mariana.gonzales/storage/outfile.csv")