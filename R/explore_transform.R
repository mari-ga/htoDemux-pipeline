#!/usr/bin/env Rscript
#Receive arguments from command line

options(echo=TRUE)


#Libraries
library(Seurat)
library(DropletUtils, warn.conflicts = FALSE)
library(argparser, quietly=TRUE)

p <- arg_parser("Parameters for HTODemux")

p <- add_argument(p, "--fileUmi",help="Path to file UMI count matrix")
p <- add_argument(p, "--fileHto",help="Path to file HTO matrix")
p <- add_argument(p, "--path10x",help="Path to save 10x matrix")



argv <- parse_args(p)


rna <- Read10X(data.dir = argv$fileUmi)
hto <- Read10X(data.dir = argv$fileHto)
#rna <- readRDS(argv$fileUmi)
#hto <- readRDS(argv$fileHto)

typeof(hto)
print("--------------------")
dim(rna)
print("--------------------")
dim(hto)
#tmpdir <- tempfile()
#print(tmpdir)
#write10xCounts("/Users/mylenemarianagonzalesandre/Development/Data-cluster/tutorial/hto/",hto)
#transposed <-t(hto)

#create_files <- function(name,extension) {
#  path_complete <- paste(name,extension,sep="")
#  print(path_complete)
#  if (file.exists(path_complete)) {
#    print("The file already exists...")
#    return(-1)
#  } else {
#    print("Created new file with results")
#    file.create(path_complete)
#    return(path_complete)
#  }
#}

#file_results <-create_files("hto_mat",".csv")
#write.csv(transposed, file=file_results)