#!/usr/bin/env Rscript
#BiocManager::install("DropletUtils")


#Receive arguments from command line
options(echo=TRUE)

#Libraries
library(Seurat)
library(DropletUtils, warn.conflicts = FALSE)
library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("Parameters for Hashed Drops Demultiplexing")

#Parameters - section 1
#Import files

#p <- add_argument(p, "--fileUmi",help="Path to file UMI count matrix")

p <- add_argument(p, "--fileHto",help="Path to file HTO matrix")
p <- add_argument(p, "--rdsObject",help="True if inputs are rds objects")
p <- add_argument(p, "--rawData",help="True if inputs comes from raw data - previously treated with empty drops",default=FALSE)

#Output paths
#p <- add_argument(p, "--hashedDropsPath",help="Path to file where the results of Hashed Drops will be saved", default = NULL)
p <- add_argument(p, "--nameOutputFileDrops",help="Name for the file containing the output of Hashed Drop object", default = "resultHashed.csv")
p <- add_argument(p, "--nameOutputFileHashed",help="Name for the rds Object containing the Hashed Drop results", default = "resultHashed_object.rds")

#Demultiplexing parameters
p <- add_argument(p, "--ambient",help="Specifies the relative abundance of each HTO in the ambient solution", default = "NULL")
p <- add_argument(p, "--minProp",help="infer the ambient profile when ambient=NULL", default = 0.05)
p <- add_argument(p, "--pseudoCount",help="minimum pseudo-count when computing log-fold changes", default = 5)
p <- add_argument(p, "--constAmbient",help=" indicates whether a constant level of ambient contamination should be used to estimate LogFC2 for all cells", default = FALSE)
p <- add_argument(p, "--doubletNmads",help="Specifies the number of median absolute deviations (MADs) to use to identify doublets.", default = 3)
p <- add_argument(p, "--doubletMin",help="Specifies the number of median absolute deviations (MADs) to use to identify doublets.", default = 2)
p <- add_argument(p, "--confidenMin",help="Specifies the minimum threshold on the log-fold change to use to identify singlets.", default = 2)
p <- add_argument(p, "--confidentNmads",help="Specifies the number of MADs to use to identify confidently assigned singlet", default = 3)

argv <- parse_args(p)


#---------------- Section 1 - Input files -----------------

#If the data is raw, we receive the dataframe from empty drops
print(argv$rdsObject)


counts <- Read10X(data.dir = argv$fileHto)
if(isTRUE(argv$rawData)){
 
 empty_results <- readRDS(argv$rdsObject)
 #emptyDrops(counts)
 print(empty_results)
 #Code from BioConductor: https://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html#demultiplexing-hashed-libraries
 has.cell <- empty_results$FDR <= 0.001
 print(summary(has.cell))
 print(counts[,which(has.cell)])
 hashed <- hashedDrops(counts[,which(has.cell)], ambient=metadata(empty_results)$ambient,min.prop = argv$minProp, constant.ambient = argv$constAmbient,doublet.nmads=argv$doubletNmads,confident.min=argv$confidenMin,confident.nmads=argv$confidentNmads,doublet.min=argv$doubletMin)

}else{
print("Hashed Drops no raw data")
hashed <- hashedDrops(counts, ambient = NULL, min.prop = argv$minProp, constant.ambient = argv$constAmbient,doublet.nmads=argv$doubletNmads,confident.min=argv$confidenMin,confident.nmads=argv$confidentNmads,doublet.min=argv$doubletMin)
}




#------------------Section 3 - Saving results ---------------------------"

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

file_results <-create_files(argv$nameOutputFileDrops,".csv")

write.csv(hashed, file=file_results)
pbmc_file = paste(argv$nameOutputFileHashed,".rds",sep="")
saveRDS(hashed, file=pbmc_file)

