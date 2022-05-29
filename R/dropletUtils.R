#BiocManager::install("DropletUtils")


#Receive arguments from command line
options(echo=TRUE)

#Libraries
library(Seurat)
library(DropletUtils)
library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("Parameters for Seurat object")

#Parameters - section 1
#Import files

p <- add_argument(p, "--fileUmi",help="Path to file UMI count matrix")
p <- add_argument(p, "--fileHto",help="Path to file HTO matrix")

argv <- parse_args(p)


#---------------- Section 1 - Input files -----------------
pbmc.umis <-readRDS(argv$fileUmi)
pbmc.htos <-readRDS(argv$fileHto)

#Identify which UMI corresponds to which hashtag.
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Confirm that the HTO have the correct names
rownames(pbmc.htos)

hashedDrops(pbmc.htos,  ambient = NULL,min.prop = 0.05,constant.ambient = FALSE,)