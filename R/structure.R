#!/usr/bin/env Rscript
library(Seurat)
library(argparser, quietly=TRUE)
library(ggplot2)

# Create a parser
p <- arg_parser("Parameters for Seurat object")

#Parameters - section 1
#Import files

p <- add_argument(p, "--fileUmi",help="Path to file UMI count matrix")
p <- add_argument(p, "--fileHto",help="Path to file HTO matrix")
argv <- parse_args(p)


#pbmc.umis <-readRDS(argv$fileUmi)
pbmc.htos <-readRDS(argv$fileHto)

#pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis)

print(pbmc.htos)

str(pbmc.htos)