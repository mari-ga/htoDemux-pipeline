#!/usr/bin/env Rscript
#install.packages("Seurat",repos=("http://cran.rstudio.com"))
#install.packages("spatstat.sparse",repos=("http://cran.rstudio.com"))
#install.packages('argparser', repos=("http://cran.rstudio.com"))
#install.packages("devtools", repos=("http://cran.rstudio.com"))
#install.packages("data.table", repos=("http://cran.rstudio.com"))


#Receive arguments from command line
options(echo=TRUE)

#Libraries
library(Seurat)
library(argparser, quietly=TRUE)



# Create a parser
p <- arg_parser("Parameters for Seurat object")

#Parameters - section 1
#Import files

p <- add_argument(p, "--fileUmi",help="Path to file UMI count matrix")
p <- add_argument(p, "--fileHto",help="Path to file HTO matrix")
p <- add_argument(p, "--ndelim",help="For the initial identity calss for each cell, delimiter for the cell's column name",default="_")
p <- add_argument(p, "--rdsObject",help="True if inputs are rds objects",default=FALSE)

#Parameters - section 2
p <- add_argument(p, "--selectMethod",help="Selection method", default="mean.var.plot")
p <- add_argument(p, "--numberFeatures",help="Number of features to be used when finding variable features", type="numeric", default=2000)
p <- add_argument(p, "--assay",help="Choose assay between RNA or HTO",default="HTO")

#Parameters - section 3
p <- add_argument(p, "--normalisationMethod",help="Normalisation method", default="CLR")
p <- add_argument(p, "--margin",help="Margin for normalisation", type="numeric",default=2)
p <- add_argument(p, "--assayName",help="Name of the Hashtag assay HTO by default",default="HTO")

#Output paths

p <- add_argument(p, "--nameOutputFile",help="Name for the file containing the output of HTODemux hashtag", default = "object")

argv <- parse_args(p)

#---------------- Section 1 - Input files -----------------
#umi stands for the RNA matrix
if(isTRUE(argv$rdsObject)){
umi <- readRDS(argv$fileUmi)
counts <- readRDS(argv$fileHto)
}else{
umi <- Read10X(data.dir = argv$fileUmi)
counts <- Read10X(data.dir = argv$fileHto)
}



#Identify which UMI corresponds to which hashtag.
joint.bcs <- intersect(colnames(umi), colnames(counts))
print(joint.bcs)


umi<- umi[, joint.bcs]
counts <- counts[, joint.bcs]

#print(joint.bcs)
# Subset RNA and HTO counts by joint cell barcodes
#pbmc.umis <- umi[, joint.bcs]
#pbmc.htos <- as.matrix(counts[, joint.bcs])

# Confirm that the HTO have the correct names
#rownames(pbmc.htos)

#-------------------- Section 2 - Setup Seurat ---------------------------------------

#Setup Seurat object and add in the HTO data
# Setup Seurat object
pbmc.hashtag <- CreateSeuratObject(counts = umi, names.delim = argv$ndelim)
#For RNA data
#pbmc.hashtag <-FindVariableFeatures(pbmc.hashtag, selection.method=argv$selectMethod)
#pbmc.hashtag <- ScaleData(pbmc.hashtag, features=  VariableFeatures(pbmc.hashtag))

#------------------ Section 3 - adding HTO data as an independent assay ---------------------
# Add HTO data as a new assay independent from RNA
pbmc.hashtag[[argv$assayName]] <- CreateAssayObject(counts = counts)
# Normalize HTO data
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = argv$assayName, normalization.method = argv$normalisationMethod, margin=argv$margin)


#------------------Section 5 - Save object for demultiplex ---------------------------

create_files <- function(name, path,extension,converter) {
  path_complete <- paste(path, name,extension,sep="")
  print(path_complete)
  if(isTRUE(converter))
  {
    if (file.exists(path_complete)) {
      print("The file already exists...")
      return(-1)
    } else {
      print("Created new file with results")
      file.create(path_complete)
      return(path_complete)
    }
  }else{
    print("We don't need to convert the input files")
  }
}


#Save Results
#pbmc_file = paste(argv$demulOutPath,argv$nameOutputFile,".rds",sep="")
pbmc_file = paste(argv$nameOutputFile,".rds",sep="")
saveRDS(pbmc.hashtag, file=pbmc_file)





