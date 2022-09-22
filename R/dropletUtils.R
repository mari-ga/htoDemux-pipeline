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

p <- add_argument(p, "--fileUmi",help="Path to file UMI count matrix")
p <- add_argument(p, "--fileHto",help="Path to file HTO matrix")


#Output paths
#p <- add_argument(p, "--hashedDropsPath",help="Path to file where the results of Hashed Drops will be saved", default = NULL)
p <- add_argument(p, "--nameOutputFileDrops",help="Name for the file containing the output of Hashed Drop object", default = "resultHashed.csv")
p <- add_argument(p, "--nameOutputFileHashed",help="Name for the rds Object containing the Hashed Drop results", default = "resultHashed.rds")

#Demultiplexing parameters
p <- add_argument(p, "--ambient",help="Specifies the relative abundance of each HTO in the ambient solution", default = "NULL")
p <- add_argument(p, "--minProp",help="infer the ambient profile when ambient=NULL", default = 0.05)
p <- add_argument(p, "--pseudoCount",help="minimum pseudo-count when computing log-fold changes", default = 5)
p <- add_argument(p, "--constAmbient",help=" indicates whether a constant level of ambient contamination should be used to estimate LogFC2 for all cells", default = FALSE)
p <- add_argument(p, "--doubletNmads",help="Specifies the number of median absolute deviations (MADs) to use to identify doublets.", default = 3)
p <- add_argument(p, "--doubletMin",help="Specifies the number of median absolute deviations (MADs) to use to identify doublets.", default = 2)
p <- add_argument(p, "--confidenMin",help="Specifies the minimum threshold on the log-fold change to use to identify singlets.", default = 2)
p <- add_argument(p, "--confidentNmads",help="Specifies the number of MADs to use to identify confidently assigned singlet", default = 3)

p <- add_argument(p, "--combinations",help="Specifies valid combinations of HTOs", default = NULL)

argv <- parse_args(p)


print(argv$ambient)
#---------------- Section 1 - Input files -----------------
# pbmc.umis <-readRDS(argv$fileUmi)
# pbmc.htos <-readRDS(argv$fileHto)

# #Identify which UMI corresponds to which hashtag.
# joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# # Subset RNA and HTO counts by joint cell barcodes
# pbmc.umis <- pbmc.umis[, joint.bcs]
# pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# # Confirm that the HTO have the correct names
# rownames(pbmc.htos)
# print(argv$ambient)


umi <- Read10X(data.dir = argv$fileUmi)
counts <- Read10X(data.dir = argv$fileHto)


#Identify which UMI corresponds to which hashtag.
joint.bcs <- intersect(colnames(umi), colnames(counts))

# Subset RNA and HTO counts by joint cell barcodes
umi<- umi[, joint.bcs]
counts <- counts[, joint.bcs]

#---------------- Section 2 - Demultiplexing -----------------
#hashed <- hashedDrops(pbmc.htos,  ambient = argv$ambient ,min.prop = argv$minProp, pseudo.count=argv$pseudoCount, constant.ambient = argv$constAmbient, doublet.nmads=argv$doubletNmads, confident.min=argv$confidenMin ,combinations=argv$combinations,confident.nmads=argv$confidentNmads,doublet.min=arg$doubletMin)
#hashed <- hashedDrops(pbmc.htos,  ambient = argv$ambient ,min.prop = argv$minProp, constant.ambient = argv$constAmbient)
hashed <- hashedDrops(counts, ambient = NULL, min.prop = argv$minProp, constant.ambient = argv$constAmbient,doublet.nmads=argv$doubletNmads)

#hashed <- hashedDrops(pbmc.umis,  ambient = NULL,min.prop = 0.05,constant.ambient = FALSE,)

str(hashed)

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
print(pbmc_file)
saveRDS(hashed, file=pbmc_file)

