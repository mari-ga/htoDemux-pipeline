#!/usr/bin/env Rscript
#install.packages("Seurat",repos=("http://cran.rstudio.com"))
#install.packages("spatstat.sparse",repos=("http://cran.rstudio.com"))
#install.packages('argparser', repos=("http://cran.rstudio.com"))
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

p <- add_argument(p, "--fileUmi",help="Path to file UMI count matrix")
p <- add_argument(p, "--fileHto",help="Path to file HTO matrix")

#Parameters - section 2
p <- add_argument(p, "--selectMethod",help="Selection method", default="vst")
p <- add_argument(p, "--numberFeatures",help="Number of features to be used when finding variable features", type="numeric", default=2000)
p <- add_argument(p, "--assay",help="Choose assay between RNA or HTO",default="HTO")

#Parameters - section 3
p <- add_argument(p, "--normalisationMethod",help="Normalisation method", default="CLR")
p <- add_argument(p, "--margin",help="Margin for normalisation", type="numeric",default=2)
p <- add_argument(p, "--assayName",help="Name of the Hashtag assay HTO by default",default="HTO")

#parameters - section 4
p <- add_argument(p, "--quantile",help="Positive quantile per default: 0.99", type="numeric",default=0.99)
p <- add_argument(p, "--kfunc",help="Cluster function choose between: Clara - kmeans",default="clara")
p <- add_argument(p, "--nstarts",help="number of starts for demultiplex", type="numeric",default=100)
p <- add_argument(p, "--nsamples",help="number of samples for demultiplex", type="numeric",default=100)

#Output paths
p <- add_argument(p, "--htoDemuxOut",help="Path to file where the results of htoDemux will be saved", default = NULL)
p <- add_argument(p, "--nameOutputFile",help="Name for the file containing the output of HTODemux hashtag", default = "result.csv")
p <- add_argument(p, "--graphs",help="Path to folder where the graphs produced from the HTODemux function can be saved", default = NULL)
#Output graphs - Ridge Plot
p <- add_argument(p, "--ridgePlot",help="Generates a ridge plot from the results, True to generate", default = FALSE)
p <- add_argument(p, "--ridgeNCol",help="Number of columns for ridgePlot", default = 2)
#Output graphs - Scatter Feature
p <- add_argument(p, "--featureScatter",help="Generates a ridge plot from the results, True to generate", default = FALSE)
p <- add_argument(p, "--scatterFeat1",help="Feature 1 for Feature Scatter Plot", default = NULL)
p <- add_argument(p, "--scatterFeat2",help="Feature 2 for Feature Scatter Plot", default = NULL)
#Output graphs - Violin Plot
p <- add_argument(p, "--vlnplot",help="Generates a violin plot from the results, True to generate", default = FALSE)
p <- add_argument(p, "--vlnFeatures",help="Features to plot (gene expression, metrics, PC scores, anything that can be retreived by FetchData)", default = NULL)
p <- add_argument(p, "--vlnLog",help="plot the feature axis on log scale", default = TRUE)

#Output graphs - tSNE
p <- add_argument(p, "--tsne",help="Generate a two dimensional tSNE embedding for HTOs", default = FALSE)
p <- add_argument(p, "--tseIdents",help="What should we remove from the object (we have Singlet,Doublet and Negative)", default = "Negative")
p <- add_argument(p, "--tsneInvert",help="True or False", default = TRUE)
p <- add_argument(p, "--tsePerplexity",help="value for perplexity", type="numeric",  default = 100)

#Output graphs - Heatmap
p <- add_argument(p, "--heatmap",help="Generate a Heatmap", default = FALSE)
p <- add_argument(p, "--heatmapNcells",help="value for number of cells", type="numeric",  default = 5000)

#Output graphs - Cluster and PCA
p <- add_argument(p, "--cluster",help="Cluster and visualize cells - perform PCA", default = FALSE)
p <- add_argument(p, "--clusterIdents",help="Could choose Singlet Doublet", default = "Singlet")
p <- add_argument(p, "--clusterSelMethod",help="Selection method for cluster, same options available as selection method for main experiment", default = "mean.var.plot")
p <- add_argument(p, "--reductionMethod",help="Methos used to reduce the data", default = "pca")
p <- add_argument(p, "--reductionDims",help="Dimensions of reduction to use as input", type="numeric",  default = 2)
p <- add_argument(p, "--reductionResol",help="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities", type="numeric",  default = 0.1)
p <- add_argument(p, "--dimPlot",help="Generate Dimensional reduction plot", default = FALSE)


argv <- parse_args(p)


#---------------- Section 1 - Input files -----------------
pbmc.umis <-readRDS(argv$fileUmi)
print(pbmc.umis)

pbmc.htos <- readRDS(argv$fileHto)
#pbmc.htos

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
pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis, assay =argv$assay )
pbmc.hashtag

# Normalize RNA data with log normalization
pbmc.hashtag <- NormalizeData(pbmc.hashtag)
# Find and scale variable features
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = argv$selectMethod, nfeatures=argv$numberFeatures)
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))

#------------------ Section 3 - adding HTO data as an independent assay ---------------------
# Add HTO data as a new assay independent from RNA
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = argv$assayName, normalization.method = argv$normalisationMethod, margin=argv$margin)


#------------------ Section 4 - Demultiplex cells based on HTO enrichment ---------------------

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = argv$assayName, positive.quantile = argv$quantile, kfunc = argv$kfunc, nstarts = argv$nstarts, nsamples = argv$nsamples)


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



file_results <-create_files(argv$nameOutputFile, argv$htoDemuxOut,".csv")
write.csv(pbmc.hashtag$HTO_classification.global, file=file_results)




#------------------ Section 6 - Visualisation ---------------------
graphsPath <- argv$graphs
# Ridge Plot
# Group cells based on the max HTO signal
if(argv$ridgePlot){
  Idents(pbmc.hashtag) <- "HTO_maxID"
  plot<-RidgePlot(pbmc.hashtag, assay = argv$assayName, features = rownames(pbmc.hashtag[[argv$assayName]])[1:2], ncol = argv$ridgeNCol)
  png(paste(graphsPath,"ridge.png",sep=""))
  print(plot)
  dev.off()
}

if(argv$featureScatter){
  plot2<- FeatureScatter(pbmc.hashtag, feature1 = argv$scatterFeat1 , feature2 = argv$scatterFeat2)
  png(paste(graphsPath,"FeatureScatter.png",sep=""))
  print(plot2)
  dev.off()
}

if(argv$vlnplot){
  Idents(pbmc.hashtag) <- "HTO_classification.global"
  plot3<-VlnPlot(pbmc.hashtag, features = pbmc.hashtag$nCount_RNA, pt.size = 0.1, log = TRUE)
  png(paste(graphsPath,"violinPlot.png",sep=""))
  print(plot3)
  dev.off()
}



if(argv$tsne){
  # First, we will remove negative cells from the object
  pbmc.hashtag.subset <- subset(pbmc.hashtag, idents = "Negative", invert = TRUE)
  
  # Calculate a distance matrix using HTO
  hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = pbmc.hashtag.subset, assay = "HTO"))))
  
  # Calculate tSNE embeddings with a distance matrix
  pbmc.hashtag.subset <- RunTSNE(pbmc.hashtag.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
  plot4<-DimPlot(pbmc.hashtag.subset)
  png(paste(graphsPath,"tSNE.png",sep=""))
  print(plot4)
  dev.off()
}

if(argv$heatmap){
  # To increase the efficiency of plotting, you can subsample cells using the num.cells argument
  plot5 <-HTOHeatmap(pbmc.hashtag, assay = argv$assayName, ncells = argv$heatmapNcells)
  png(paste(graphsPath,"heatmap.png",sep=""))
  print(plot5)
  dev.off()
}

if(argv$cluster){
  pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")
  pbmc.singlet <- FindVariableFeatures(pbmc.singlet, selection.method = argv$clusterSelMethod )
  pbmc.singlet <- ScaleData(pbmc.singlet, features = VariableFeatures(pbmc.singlet))
  pbmc.singlet <- RunPCA(pbmc.singlet, features = VariableFeatures(pbmc.singlet))
  
  pbmc.singlet <- FindNeighbors(pbmc.singlet, reduction = argv$reductionMethod, dims = 1:argv$reductionDims)
  pbmc.singlet <- FindClusters(pbmc.singlet, resolution = argv$reductionResol)
  
  pbmc.singlet <- RunTSNE(pbmc.singlet, reduction = argv$reductionMethod, dims = 1:argv$reductionDims)
  if(argv$dimPlot){
    plot6 <-DimPlot(pbmc.singlet, group.by = "HTO_classification")
    png(paste(graphsPath,"dimPlot.png",sep=""))
    print(plot6)
    dev.off()
  }
  
  
}



