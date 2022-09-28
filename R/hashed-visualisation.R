#!/usr/bin/env Rscript
#BiocManager::install("DropletUtils")


#Receive arguments from command line
options(echo=TRUE)

#Libraries
library(Seurat)
library(DropletUtils, warn.conflicts = FALSE)
library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("Parameters for Hashed Drops - Visualisation")

p <- add_argument(p, "--hashedObjectPath",help="hashed drops object containing results", default = NULL)
p <- add_argument(p, "--histogram",help="Plot histogram for Hashed Log2", default = TRUE)
p <- add_argument(p, "--plotLog",help="Plot Log 2 FC", default = TRUE)

argv <- parse_args(p)

hashed <-readRDS(argv$hashedObjectPath)
str(hashed)
# Doublets show up in the top-left, singlets in the bottom right.
if(argv$histogram){
plot1 <- plot(hashed$LogFC, hashed$LogFC2)
png(paste("plot_hashed.png",sep=""))
print(plot1)
}

# Most cells should be singlets with low second log-fold changes.
if(argv$plotLog){
plot2 <- hist(hashed$LogFC2, breaks=50)
png(paste("histo_hashed.png",sep=""))
print(plot2)
}
# Identify confident singlets or doublets at the given threshold.
#summary(hashed$Confident)
#summary(hashed$Doublet)

# Checking against the known truth, in this case
# 'Best' contains the putative sample of origin.
#table(hashed$Best, true.sample) 
