#!/usr/bin/env Rscript
#install.packages("Seurat",repos=("http://cran.rstudio.com"))
#install.packages("spatstat.sparse",repos=("http://cran.rstudio.com"))
#Receive arguments from command line
options(echo=TRUE)
#Get arguments as a vector
myargs = commandArgs(trailingOnly=TRUE)

#Library
library(Seurat)

#Import files
file_umis = myargs[1]
file_htos = myargs[2]


#Parameters - section 2
selection_method = myargs[3]
number_features = myargs[4]
typeof(number_features)


if(is.na(selection_method))
{
  print("Empty args for selection method")
  selection_method = "vst"
}
if(is.na(number_features))
{
  print("Empty args for number of features")
  number_features = 2000
}

#Parameters - section 3
normalisation_method = myargs[5]
margin  = myargs[6]
assay   = myargs[7]
print(assay)
if(is.na(normalisation_method))
{
  print("Empty args for normalisation method")
  print("Normalisation method per default CLR")
  normalisation_method = "CLR"
}

if(is.na(margin))
{
  print("Empty args for margin")
  margin = 2
}

if(is.na(assay))
{
  print("empty args for assay")
  print("Assay per default: HTO")
  assay  = "HTO"
}

#parameters - section 4
positive_quantile = myargs[8]
kfunc = myargs[9]
nstarts = myargs[10]
nsamples = myargs[11]

if(is.na(positive_quantile))
{
  print("Positive quantile per default: 0.99")
  positive_quantile  = 0.99
}

if(is.na(kfunc))
{
  print("Clustering method per default: Kmeans")
  kfunc  = "kmeans"
}

if(is.na(nstarts))
{
  print("Empty args for nstars, value per default")
  nstarts = 100
}

if(is.na(nsamples))
{
  nsamples = 100
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
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))

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
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = assay, normalization.method = normalisation_method, margin=margin)

#------------------ Section 4 - Demultiplex cells based on HTO enrichment ---------------------

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = assay, positive.quantile = positive_quantile, kfunc = kfunc, nstarts = nstarts, nsamples = nsamples)

# Global classification results
table(pbmc.hashtag$HTO_classification.global)

