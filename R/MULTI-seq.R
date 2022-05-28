#!/usr/bin/env Rscript
#install.packages("Seurat",repos=("http://cran.rstudio.com"))
#install.packages("spatstat.sparse",repos=("http://cran.rstudio.com"))
#install.packages('argparser', repos=("http://cran.rstudio.com"))
#Receive arguments from command line
options(echo=TRUE)

#Libraries
library(Seurat)
library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("Parameters for MULTI-seq")

p <- add_argument(p, "--seuratObjectPath",help="seurat object ready for demultiplex step", default = NULL)



#parameters - section 4
p <- add_argument(p, "--quantile",help="Positive quantile per default: 0.7", type="numeric",default=0.7)
p <- add_argument(p, "--autoThresh",help="Whether to perform automated threshold finding to define the best quantile, per Default=False",default=TRUE)
p <- add_argument(p, "--maxiter",help="Maximum number of iterations", type="numeric",default=5)
p <- add_argument(p, "--qrangeFrom",help="A range of possible quantile values to try if autoThresh is TRUE",type="numeric", default = 0.1)
p <- add_argument(p, "--qrangeTo",help="A range of possible quantile values to try if autoThresh is TRUE",type="numeric", default = 0.9)
p <- add_argument(p, "--qrangeBy",help="A range of possible quantile values to try if autoThresh is TRUE",type="numeric", default = 0.05)
p <- add_argument(p, "--verbose",help="Prints the output", default = TRUE)

#Output paths
p <- add_argument(p, "--multiSeqOutPath",help="Path to file where the results of MULTI-seq will be saved", default = NULL)
p <- add_argument(p, "--nameOutputFileMulti",help="Name for the file containing the output of MULTI-Seq object", default = "result.csv")

argv <- parse_args(p)



pbmc.hashtag <-readRDS(argv$seuratObjectPath)
str(pbmc.hashtag)


#------------------ Section 4 - Demultiplex cells based on HTO enrichment ---------------------

pbmc.hashtag <- MULTIseqDemux(pbmc.hashtag, assay = argv$assayName,  quantile = args$quantile, autoThresh = TRUE , qrange=seq(from = argv$qrangeFrom, to =argv$qrangeTo, by=argv$qrangeBy), verbose=argv$verbose)


table(pbmc.hashtag$MULTI_ID)
print("----------------------------------------------------------------------------")
table(pbmc.hashtag$MULTI_classification)
print("----------------------------------------------------------------------------")


print("-----------------------------------------------")

pbmc.hashtag
print("-----------------------------------------------")

dim(x = pbmc.hashtag)
print("-----------------------------------------------")

head(x = rownames(x = pbmc.hashtag))

head(x = colnames(x = pbmc.hashtag))

names(x = pbmc.hashtag)

print("-----------------------------------------------")

pbmc.hashtag[['RNA']]
print("-----------------------------------------------")
pbmc.hashtag[['HTO']]

print("-----------------------------------------------")

colnames(x = pbmc.hashtag[[]])


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
print(argv$nameOutputFileMulti)
print("-------")
file_results <-create_files(argv$nameOutputFileMulti, argv$multiSeqOutPath,".csv")
write.csv(pbmc.hashtag$MULTI_ID, file=file_results)
pbmc_file = paste(argv$multiSeqOutPath,argv$nameOutputFileMulti,".rds",sep="")
print(pbmc_file)
saveRDS(pbmc.hashtag, file=pbmc_file)
