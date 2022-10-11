#!/usr/bin/env Rscript
#Receive arguments from command line
options(echo=TRUE)

#Libraries
library(Seurat)
library(DropletUtils, warn.conflicts = FALSE)
library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("Parameters for Empty Drops cell identification")

p <- add_argument(p, "--fileUmi",help="Path to file UMI count matrix")

p <- add_argument(p, "--niters",help="An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations.",default=10000)
p <- add_argument(p, "--empty",help="True only if the data provided is RAW", default = TRUE)
p <- add_argument(p, "--lower",help=" numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.", default = 100)
p <- add_argument(p, "--testAmbient",help=" logical scalar indicating whether results should be returned for barcodes with totals less than or equal to lower", default = FALSE)
p <- add_argument(p, "--alpha",help="A numeric scalar specifying the scaling parameter for the Dirichlet-multinomial sampling scheme", default = NULL)
p <- add_argument(p, "--ignore",help="A numeric scalar specifying the lower bound on the total UMI count, at or below which barcodes will be ignored", default = NULL)

p <- add_argument(p, "--nameOutputEmpty",help="Name for the empty droplets file", default = "emptyDropletsHashed")

argv <- parse_args(p)

umi <- Read10X(data.dir = argv$fileUmi)
print(argv$alpha)
if(isTRUE(argv$empty)){
emptyHashed <- testEmptyDrops(umi,lower= argv$lower,niters =argv$niters,test.ambient =argv$testAmbient )
#alpha produces an error when using default values
#alpha=argv$alpha,
#Currently ignore is not supported
#ignore =argv$ignore
}
create_files <- function(name,extension) {
  path_complete <- paste( name,extension,sep="")
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

if(isTRUE(argv$empty))
{
  empty_results <-create_files(argv$nameOutputEmpty,".csv")
  write.csv(emptyHashed, file=empty_results)
}