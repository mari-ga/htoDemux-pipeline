#!/usr/bin/env Rscript
#install.packages("Seurat",repos=("http://cran.rstudio.com"))
#install.packages("spatstat.sparse",repos=("http://cran.rstudio.com"))
#install.packages('argparser', repos=("http://cran.rstudio.com"))
#install.packages("devtools", repos=("http://cran.rstudio.com"))
#install.packages("BiocManager",repos=("http://cran.rstudio.com"))
#BiocManager::install("rhdf5")
#devtools::install_github("AllenInstitute/scrattch.io")


#Receive arguments from command line
options(echo=TRUE)

#Libraries
library(Seurat)
library(argparser, quietly=TRUE)
library(ggplot2)
library(scrattch.io)
library(rhdf5)

# Create a parser
p <- arg_parser("Parameters for Seurat object")

#Parameters - section 1
#Import files

p <- add_argument(p, "--fileUmi",help="Path to file UMI count matrix")
p <- add_argument(p, "--fileHto",help="Path to file HTO matrix")

#Converter for DemuxEM
p <- add_argument(p, "--converter",help="Transform both input matrices into csv and hdf5 respectively for demuxEM", default = FALSE)
p <- add_argument(p, "--conversionPath",help="Path to save both converted files", default = NULL)
#Parameters - section 2
p <- add_argument(p, "--selectMethod",help="Selection method", default="vst")
p <- add_argument(p, "--numberFeatures",help="Number of features to be used when finding variable features", type="numeric", default=2000)
p <- add_argument(p, "--assay",help="Choose assay between RNA or HTO",default="HTO")

#Parameters - section 3
p <- add_argument(p, "--normalisationMethod",help="Normalisation method", default="CLR")
p <- add_argument(p, "--margin",help="Margin for normalisation", type="numeric",default=2)
p <- add_argument(p, "--assayName",help="Name of the Hashtag assay HTO by default",default="HTO")

#Output paths
p <- add_argument(p, "--demulOutPath",help="Path to file where the rds object ready for demultiplexing will be saved", default = NULL)
p <- add_argument(p, "--nameOutputFile",help="Name for the file containing the output of HTODemux hashtag", default = "result")

argv <- parse_args(p)

#---------------- Section 1 - Input files -----------------
pbmc.umis <-readRDS(argv$fileUmi)
pbmc.htos <-readRDS(argv$fileHto)

# Unchanged original files for conversion purposes
hto <-readRDS(argv$fileHto)
umi <-readRDS(argv$fileUmi)

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

str(pbmc.hashtag)


#------------------ Section 3 - adding HTO data as an independent assay ---------------------
# Add HTO data as a new assay independent from RNA
pbmc.hashtag[[argv$assayName]] <- CreateAssayObject(counts = pbmc.htos)
# Normalize HTO data
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = argv$assayName, normalization.method = argv$normalisationMethod, margin=argv$margin)

str(pbmc.hashtag)

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
print(argv$nameOutputFile)
pbmc_file = paste(argv$demulOutPath,argv$nameOutputFile,".rds",sep="")
print(pbmc_file)
saveRDS(pbmc.hashtag, file=pbmc_file)

#-------------- Section 6 - convert files for demuxEM (optional) --------------------------

hto_file = create_files("hto_matrix", argv$conversionPath,".csv",argv$converter)
if (hto_file != -1){
  write.csv(hto, file=hto_file) 
}else{
  print("It wasn't possible to create csv from HTO matrix")
}

if (isTRUE(argv$converter)){
  umi_file = paste(argv$conversionPath,"umi_matrix.h5",sep="")
  #h5createFile(umi_file)
  write_dgCMatrix_h5(umi, cols_are = "gene_names", umi_file,
                     ref_name = "conversion?", gene_ids = NULL)
}else{
  print("It wasn't possible to create h5 from umi matrix")
}



