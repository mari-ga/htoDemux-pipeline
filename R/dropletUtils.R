#BiocManager::install("DropletUtils")


#Receive arguments from command line
options(echo=TRUE)

#Libraries
library(Seurat)
library(DropletUtils)
library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("Parameters for Seurat object")

#Parameters - section 1
#Import files

p <- add_argument(p, "--fileUmi",help="Path to file UMI count matrix")
p <- add_argument(p, "--fileHto",help="Path to file HTO matrix")


#Output paths
p <- add_argument(p, "--hashedDropsPath",help="Path to file where the results of Hashed Drops will be saved", default = NULL)
p <- add_argument(p, "--nameOutputFileDrops",help="Name for the file containing the output of Hashed Drop object", default = "result.csv")

argv <- parse_args(p)


#---------------- Section 1 - Input files -----------------
pbmc.umis <-readRDS(argv$fileUmi)
pbmc.htos <-readRDS(argv$fileHto)

#Identify which UMI corresponds to which hashtag.
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Confirm that the HTO have the correct names
rownames(pbmc.htos)

#---------------- Section 2 - Demultiplexing -----------------
hashed <- hashedDrops(pbmc.htos,  ambient = NULL,min.prop = 0.05,constant.ambient = FALSE,)
#hashed <- hashedDrops(pbmc.umis,  ambient = NULL,min.prop = 0.05,constant.ambient = FALSE,)
str(hashed)

print("------------------")
typeof(hashed)

print("------------------")
str(pbmc.htos)



#------------------Section 3 - Saving results ---------------------------"

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

file_results <-create_files(argv$nameOutputFileDrops, argv$hashedDropsPath,".csv")
write.csv(hashed, file=file_results)

