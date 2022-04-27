#!/usr/bin/env Rscript
#install.packages("Seurat",repos=("http://cran.rstudio.com"))
#install.packages("spatstat.sparse",repos=("http://cran.rstudio.com"))
#install.packages('argparser', repos=("http://cran.rstudio.com"))
#Receive arguments from command line
options(echo=TRUE)

#Libraries
library(Seurat)
library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("Parameters for HTODemux")

#Parameters - section 1
#Import files

p <- add_argument(p, "--fileUmi",help="Path to file UMI count matrix")
p <- add_argument(p, "--fileHto",help="Path to file HTO matrix")

#Parameters - section 2
p <- add_argument(p, "--selectMethod",help="Selection method", default="vst")
p <- add_argument(p, "--numberFeatures",help="Number of features to be used when finding variable features", type="numeric", default=2000)

#Parameters - section 3
p <- add_argument(p, "--normalisationMethod",help="Normalisation method", default="CLR")
p <- add_argument(p, "--margin",help="Margin for normalisation", type="numeric",default=2)
p <- add_argument(p, "--assay",help="Choose assay between RNA or HTO",default="HTO")

#parameters - section 4
p <- add_argument(p, "--quantile",help="Positive quantile per default: 0.99", type="numeric",default=0.99)
p <- add_argument(p, "--kfunc",help="Cluster function choose between: Clara - kmeans",default="kmeans")
p <- add_argument(p, "--nstarts",help="number of starts for demultiplex", type="numeric",default=100)
p <- add_argument(p, "--nsamples",help="number of samples for demultiplex", type="numeric",default=100)

argv <- parse_args(p)


#---------------- Section 1 - Input files -----------------
pbmc.umis <-readRDS(argv$fileUmi)
print(pbmc.umis)

pbmc.htos <- readRDS(argv$fileHto)
#pbmc.htos

#Identify which UMI corresponds to which hashtag.
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Confirm that the HTO have the correct names
rownames(pbmc.htos)

#-------------------- Section 2 - Setup Seurat ---------------------------------------

#Setup Seurat object and add in the HTO data

# Setup Seurat object
pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis)
pbmc.hashtag

# Normalize RNA data with log normalization
pbmc.hashtag <- NormalizeData(pbmc.hashtag)
# Find and scale variable features
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = argv$selectMethod, nfeatures=argv$numberFeatures)
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))

#------------------ Section 3 - adding HTO data as an independent assay ---------------------
# Add HTO data as a new assay independent from RNA
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = argv$assay, normalization.method = argv$normalisationMethod, margin=argv$margin)


#------------------ Section 4 - Demultiplex cells based on HTO enrichment ---------------------

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = argv$assay, positive.quantile = argv$quantile, kfunc = argv$kfunc, nstarts = argv$nstarts, nsamples = argv$nsamples)

# Global classification results
table(pbmc.hashtag$HTO_classification.global)
