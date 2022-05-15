#!/usr/bin/env Rscript
install.packages("Seurat",repos=("http://cran.rstudio.com"))
install.packages("spatstat.sparse",repos=("http://cran.rstudio.com"))
install.packages('argparser', repos=("http://cran.rstudio.com"))
#Receive arguments from command line
options(echo=TRUE)

#Libraries
library(Seurat)
library(argparser, quietly=TRUE)
library(ggplot2)

# Create a parser
p <- arg_parser("Parameters for HTODemux")

#Parameters - section 1
#Import files

p <- add_argument(p, "--seuratObjectPath",help="seurat object ready for demultiplex step", default = NULL)

#parameters - section 4
p <- add_argument(p, "--quantile",help="Positive quantile per default: 0.99", type="numeric",default=0.99)
p <- add_argument(p, "--kfunc",help="Cluster function choose between: Clara - kmeans",default="clara")
p <- add_argument(p, "--nstarts",help="number of starts for demultiplex", type="numeric",default=100)
p <- add_argument(p, "--nsamples",help="number of samples for demultiplex", type="numeric",default=100)

#Output paths
p <- add_argument(p, "--htoDemuxOutPath",help="Path to file where the results of htoDemux will be saved", default = NULL)
p <- add_argument(p, "--nameOutputFile",help="Name for the file containing the output of HTODemux hashtag", default = "result.csv")


argv <- parse_args(p)

#------------------ Loading Seurat object ---------------------

pbmc.hashtag <-readRDS(argv$seuratObjectPath)
str(pbmc.hashtag)
#------------------ Section 4 - Demultiplex cells based on HTO enrichment ---------------------

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = argv$assayName, positive.quantile = argv$quantile,  nstarts = argv$nstarts, nsamples = argv$nsamples)


# Global classification results
table(pbmc.hashtag$HTO_classification.global)

#------------------Section 5 - Saving results ---------------------------"

create_files <- function(name, path,extension) {
  path_complete <- paste(path, name,extension,sep="")
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
print(argv$nameOutputFile)
print("-------")
file_results <-create_files(argv$nameOutputFile, argv$htoDemuxOutPath,".csv")
write.csv(pbmc.hashtag$HTO_classification.global, file=file_results)
pbmc_file = paste(argv$htoDemuxOutPath,argv$nameOutputFile,".rds",sep="")
print(pbmc_file)
saveRDS(pbmc.hashtag, file=pbmc_file)


