#!/usr/bin/env Rscript

options(echo=TRUE)

#Libraries
library(Seurat)
library(argparser, quietly=TRUE)
library(ggplot2)

# Create a parser
p <- arg_parser("Parameters for HTODemux")


p <- add_argument(p, "--fileUmi",help="Path to file UMI count matrix")
p <- add_argument(p, "--fileHto",help="Path to file HTO matrix")
p <- add_argument(p, "--ndelim",help="For the initial identity calss for each cell, delimiter for the cell's column name",default="_")

#Parameters - section 2
p <- add_argument(p, "--selectMethod",help="Selection method", default="mean.var.plot")
p <- add_argument(p, "--numberFeatures",help="Number of features to be used when finding variable features", type="numeric", default=2000)
p <- add_argument(p, "--assay",help="Choose assay between RNA or HTO",default="HTO")

#Parameters - section 3
p <- add_argument(p, "--normalisationMethod",help="Normalisation method", default="LogNormalize")
p <- add_argument(p, "--margin",help="Margin for normalisation", type="numeric",default=1)
p <- add_argument(p, "--assayName",help="Name of the Hashtag assay HTO by default",default="HTO")

#parameters - section 4
p <- add_argument(p, "--quantile",help="Positive quantile per default: 0.99", type="numeric",default=0.99)
p <- add_argument(p, "--kfunc",help="Cluster function choose between: Clara - kmeans",default=NULL)
p <- add_argument(p, "--nstarts",help="number of starts for demultiplex", type="numeric",default=100)
p <- add_argument(p, "--nsamples",help="number of samples for demultiplex", type="numeric",default=100)
p <- add_argument(p, "--seed",help="sets random seed", type="numeric",default=42)
p <- add_argument(p, "--init",help="Initial number of clusters for hashtags")

#Output paths
#p <- add_argument(p, "--htoDemuxOutPath",help="Path to file where the results of htoDemux will be saved", default = NULL)
p <- add_argument(p, "--nameOutputFileHTO",help="Name for the file containing the output of HTODemux object", default = "result.csv")

argv <- parse_args(p)

umi <- Read10X(data.dir = argv$fileUmi)
counts <- Read10X(data.dir = argv$fileHto)

#Identify which UMI corresponds to which hashtag.
joint.bcs <- intersect(colnames(umi), colnames(counts))


umi<- umi[, joint.bcs]
counts <- counts[, joint.bcs]

#-------------------- Section 2 - Setup Seurat ---------------------------------------

#Setup Seurat object and add in the HTO data
pbmc.hashtag <- CreateSeuratObject(counts = umi, names.delim = argv$ndelim)

#For RNA data
#pbmc.hashtag <-FindVariableFeatures(pbmc.hashtag, selection.method=argv$selectMethod)
#pbmc.hashtag <- ScaleData(pbmc.hashtag, features=  VariableFeatures(pbmc.hashtag))

#------------------ Section 3 - adding HTO data as an independent assay ---------------------
# Add HTO data as a new assay independent from RNA
pbmc.hashtag[[argv$assayName]] <- CreateAssayObject(counts = counts)
# Normalize HTO data
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = argv$assayName, normalization.method = argv$normalisationMethod, margin=argv$margin)



pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = argv$assayName, positive.quantile = argv$quantile,  nstarts = argv$nstarts, nsamples = argv$nsamples)


#------------------Section 5 - Saving results ---------------------------"

create_files <- function(name,extension) {
  path_complete <- paste(name,extension,sep="")
  print(path_complete)
  if (file.exists(path_complete)) {
    print("The file already exists...")
    return(-1)
  } else {
    print("Created new file with results")
    file.create(path_complete)
    return(path_complete)
  }
}


#Save Results
print(argv$nameOutputFileHTO)
print("-------")
file_results <-create_files(argv$nameOutputFileHTO,".csv")
write.csv(pbmc.hashtag$HTO_classification.global, file=file_results)
pbmc_file = paste(argv$nameOutputFileHTO,".rds",sep="")
saveRDS(pbmc.hashtag, file=pbmc_file)
