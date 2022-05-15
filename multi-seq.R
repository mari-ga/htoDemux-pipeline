#!/usr/bin/env Rscript
install.packages("Seurat",repos=("http://cran.rstudio.com"))
#install.packages("spatstat.sparse",repos=("http://cran.rstudio.com"))
install.packages('argparser', repos=("http://cran.rstudio.com"))
#Receive arguments from command line
options(echo=TRUE)

#Libraries
library(Seurat)
library(argparser, quietly=TRUE)
library(ggplot2)

# Create a parser
p <- arg_parser("Parameters for MULTI-seq")

p <- add_argument(p, "--fileUmi",help="Path to file UMI count matrix")
p <- add_argument(p, "--fileHto",help="Path to file HTO matrix")

#Parameters - section 2
p <- add_argument(p, "--selectMethod",help="Selection method", default="vst")
p <- add_argument(p, "--numberFeatures",help="Number of features to be used when finding variable features", type="numeric", default=2000)
p <- add_argument(p, "--assay",help="Choose assay between RNA or HTO",default="HTO")

#Parameters - section 3
p <- add_argument(p, "--normalisationMethod",help="Normalisation method", default="CLR")
p <- add_argument(p, "--margin",help="Margin for normalisation", type="numeric",default=2)
p <- add_argument(p, "--assayName",help="Name of the Hashtag assay HTO by default",default="HTO")

#parameters - section 4
p <- add_argument(p, "--quantile",help="Positive quantile per default: 0.7", type="numeric",default=0.7)
p <- add_argument(p, "--autoThresh",help="Whether to perform automated threshold finding to define the best quantile, per Default=False",default=TRUE)
p <- add_argument(p, "--maxiter",help="Maximum number of iterations", type="numeric",default=5)
p <- add_argument(p, "--qrangeFrom",help="A range of possible quantile values to try if autoThresh is TRUE",type="numeric", default = 0.1)
p <- add_argument(p, "--qrangeTo",help="A range of possible quantile values to try if autoThresh is TRUE",type="numeric", default = 0.9)
p <- add_argument(p, "--qrangeBy",help="A range of possible quantile values to try if autoThresh is TRUE",type="numeric", default = 0.05)
p <- add_argument(p, "--verbose",help="Prints the output", default = TRUE)

argv <- parse_args(p)
#---------------- Section 1 - Input files -----------------
pbmc.umis <-readRDS(argv$fileUmi)
print(pbmc.umis)

pbmc.htos <- readRDS(argv$fileHto)


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
pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis, assay =argv$assay )
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
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = argv$assayName, normalization.method = argv$normalisationMethod, margin=argv$margin)


#------------------ Section 4 - Demultiplex cells based on HTO enrichment ---------------------

pbmc.hashtag <- MULTIseqDemux(pbmc.hashtag, assay = argv$assayName,  quantile = args$quantile, autoThresh = TRUE , qrange=seq(from = argv$qrangeFrom, to =argv$qrangeTo, by=argv$qrangeBy), verbose=argv$verbose)


