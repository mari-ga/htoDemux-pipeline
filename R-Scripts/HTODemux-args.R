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

#Definir parametros por defecto - Mejores resultados por experimento de parametros

#Parameters - section 2
if(is.null(myargs[3]))
{
  selection_method = "vst"
}else{
  selection_method = myargs[3]
}

if(is.null(myargs[4]))
{
  number_features = "2000"
}else{
  number_features = myargs[4]
}
#Parameters - section 3
if(is.null(myargs[5]))
{
  print("Normalisation method per default CLR")
  normalization_method = "CLR"
}else{
  normalization_method  = myargs[5]
}

if(is.null(myargs[6]))
{
  margin = "2"
}else{
  margin  = myargs[6]
}

if(is.null(myargs[7]))
{
  print("Assay per default: HTO")
  assay  == "HTO"
}else{
  assay   = myargs[7]
}

#parameters - section 4
if(is.null(myargs[8]))
{
  print("Positive quantile per default: 0.99")
  positive_quantile  = 0.99
}else{
  positive_quantile   = myargs[8]
}
if(is.null(myargs[8]))
{
  print("Clustering method per default: Kmeans")
  kfunc  = "kmeans"
}else{
  kfunc = myargs[9]
}

if(is.null(myargs[10]))
{
  nstarts = 100
}else{
  nstarts = myargs[10]
}


if(is.null(myargs[11]))
{
  nsamples = 100
}else{
  nsamples = myargs[11]
}


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

#-------------------- Section 2 - Setup Seurat ---------------------------------------

#Setup Seurat object and add in the HTO data

# Setup Seurat object
pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis)
pbmc.hashtag

# Normalize RNA data with log normalization
pbmc.hashtag <- NormalizeData(pbmc.hashtag)
# Find and scale variable features
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = selection_method, nfeatures=number_features)
print("Using vst as selection method")
pbmc.hashtag
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))

#------------------ Section 3 - adding HTO data as an independent assay ---------------------
# Add HTO data as a new assay independent from RNA
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = assay, normalization.method = normalization_method, margin=margin)


#------------------ Section 4 - Demultiplex cells based on HTO enrichment ---------------------

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = assay, positive.quantile = positive_quantile, kfunc = kfunc, nstarts = nstarts, nsamples = nsamples)

