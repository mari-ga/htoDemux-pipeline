#!/usr/bin/env Rscript

options(echo=TRUE)

#Libraries
library(Seurat)
library(argparser, quietly=TRUE)
library(ggplot2)

# Create a parser
p <- arg_parser("Parameters for Multi-seq")


p <- add_argument(p, "--fileUmi",help="Path to file UMI count matrix")
p <- add_argument(p, "--fileHto",help="Path to file HTO matrix")
p <- add_argument(p, "--ndelim",help="For the initial identity calss for each cell, delimiter for the cell's column name",default="_")
p <- add_argument(p, "--rdsObject",help="True if inputs are rds objects",default=FALSE)

#Parameters - section 2
p <- add_argument(p, "--selectMethod",help="Selection method", default="mean.var.plot")
p <- add_argument(p, "--numberFeatures",help="Number of features to be used when finding variable features", type="numeric", default=2000)
p <- add_argument(p, "--assay",help="Choose assay between RNA or HTO",default="HTO")

#Parameters - section 3
p <- add_argument(p, "--normalisationMethod",help="Normalisation method", default="LogNormalize")
p <- add_argument(p, "--margin",help="Margin for normalisation", type="numeric",default=1)
p <- add_argument(p, "--assayName",help="Name of the Hashtag assay HTO by default",default="HTO")

#parameters - section 4
p <- add_argument(p, "--quantile_multi",help="Positive quantile per default: 0.7", type="numeric",default=0.7)
p <- add_argument(p, "--autoThresh",help="Whether to perform automated threshold finding to define the best quantile, per Default=False",default=FALSE)
p <- add_argument(p, "--maxiter",help="Maximum number of iterations", type="numeric",default=5)
p <- add_argument(p, "--qrangeFrom",help="A range of possible quantile values to try if autoThresh is TRUE",type="numeric", default = 0.1)
p <- add_argument(p, "--qrangeTo",help="A range of possible quantile values to try if autoThresh is TRUE",type="numeric", default = 0.9)
p <- add_argument(p, "--qrangeBy",help="A range of possible quantile values to try if autoThresh is TRUE",type="numeric", default = 0.05)
p <- add_argument(p, "--verbose",help="Prints the output", default = TRUE)

p <- add_argument(p, "--nameOutputFileMulti",help="Name for the file containing the output of MULTI-Seq object", default = "resultMulti_object")
p <- add_argument(p, "--nameClassificationFileMulti",help="Name for the file containing the classification of MULTI-Seq", default = "resultMulti")


argv <- parse_args(p)

#---------------- Section 1 - Input files -----------------



if(isTRUE(argv$rdsObject)){
#Read file from rds object
umi <- readRDS(argv$fileUmi)
counts <- readRDS(argv$fileHto)
}else{
#Read file from X10 files
umi <- Read10X(data.dir = argv$fileUmi)
counts <- Read10X(data.dir = argv$fileHto)
}


# Subset RNA and HTO counts by joint cell barcodes
#Identify which UMI corresponds to which hashtag.
joint.bcs <- intersect(colnames(umi), colnames(counts))


umi<- umi[, joint.bcs]
counts <-counts[, joint.bcs]

pbmc.hashtag <- CreateSeuratObject(counts = umi)
pbmc.hashtag[[argv$assayName]] <- CreateAssayObject(counts = counts)

pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = argv$assayName, normalization.method = argv$normalisationMethod)

pbmc.hashtag <- MULTIseqDemux(pbmc.hashtag, assay = argv$assayName,  quantile = args$quantile_multi, autoThresh = TRUE , qrange=seq(from = argv$qrangeFrom, to =argv$qrangeTo, by=argv$qrangeBy), verbose=argv$verbose)


#------------------Section 5 - Save results ---------------------------
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
print(argv$nameOutputFileMulti)
print("-------")
file_results <-create_files(argv$nameClassificationFileMulti,".csv")
write.csv(pbmc.hashtag$MULTI_ID, file=file_results)
pbmc_file = paste(argv$nameOutputFileMulti,".rds",sep="")
print(pbmc_file)
saveRDS(pbmc.hashtag, file=pbmc_file)
