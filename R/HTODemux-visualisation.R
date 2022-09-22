#------------- Seurat - Visualisation ------------
library(argparser, quietly=TRUE)
library(Seurat)
library(ggplot2)


# Create a parser
p <- arg_parser("Parameters for Seurat Visualisation")
p <- add_argument(p, "--pbcmHashtagPath",help="S4 object saved from the first part of HTODemux", default = NULL)
p <- add_argument(p, "--graphs",help="Path to folder where the graphs produced from the HTODemux function can be saved", default = NULL)
p <- add_argument(p, "--assayName",help="Name of the Hashtag assay HTO by default",default="HTO")

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
p <- add_argument(p, "--tsneVerbose",help="True or False", default = FALSE)
p <- add_argument(p, "--tsneApprox",help="True or False", default = FALSE)
p <- add_argument(p, "--tsneDimMax",help="max -> number of donors ",type="numeric", default=1)
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


print(argv$pbcmHashtagPath)
#"/Users/mylenemarianagonzalesandre/Development/Bachelor-Thesis/nextflow-files/htoDemux-pipeline/results/results.rds"
pbmc.hashtag <-readRDS(argv$pbcmHashtagPath)


#graphsPath <- argv$graphs
# Ridge Plot
# Group cells based on the max HTO signal
if(argv$ridgePlot){
  Idents(pbmc.hashtag) <- "HTO_maxID"
  plot<-RidgePlot(pbmc.hashtag, assay = argv$assayName, features = rownames(pbmc.hashtag[[argv$assayName]])[1:2], ncol = argv$ridgeNCol)
  png(paste("ridge.png",sep=""))
  print(plot)
  #dev.off()
}

if(argv$featureScatter){
  plot2<- FeatureScatter(pbmc.hashtag, feature1 = argv$scatterFeat1 , feature2 = argv$scatterFeat2)
  png(paste("FeatureScatter.png",sep=""))
  print(plot2)
  #dev.off()
}

if(argv$vlnplot){
  Idents(pbmc.hashtag) <- "HTO_classification.global"
  plot3<-VlnPlot(pbmc.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
  png(paste("violinPlot.png",sep=""))
  print(plot3)
  #dev.off()
}



if(argv$tsne){
  # remove negative cells from the object
  pbmc.hashtag.subset <- subset(pbmc.hashtag, idents = argv$tseIdents, invert = argv$tsneInvert)
  # Calculate a tSNE embedding of the HTO data
  DefaultAssay(pbmc.hashtag.subset) <- argv$assayName
  pbmc.hashtag.subset <- ScaleData(pbmc.hashtag.subset, features = rownames(pbmc.hashtag.subset),
                                   verbose = FALSE)
  pbmc.hashtag.subset <- RunPCA(pbmc.hashtag.subset, features = rownames(pbmc.hashtag.subset), approx = argv$tsneApprox)
  pbmc.hashtag.subset <- RunTSNE(pbmc.hashtag.subset, dims = 1:argv$tsneDimMax, perplexity = argv$tsePerplexity)
  plot4<-DimPlot(pbmc.hashtag.subset)
  png(paste("tSNE.png",sep=""))
  print(plot4)
  dev.off()
}

if(argv$heatmap){
  # To increase the efficiency of plotting, you can subsample cells using the num.cells argument
  plot5 <-HTOHeatmap(pbmc.hashtag, assay = argv$assayName, ncells = argv$heatmapNcells)
  png(paste("heatmap.png",sep=""))
  print(plot5)
  dev.off()
}

if(argv$cluster){
  pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")
  pbmc.singlet <- FindVariableFeatures(pbmc.singlet, selection.method = "mean.var.plot")
  # pbmc.singlet <- ScaleData(pbmc.singlet, features = VariableFeatures(pbmc.singlet))
  # pbmc.singlet <- RunPCA(pbmc.singlet, features = VariableFeatures(pbmc.singlet))
  # 
  # pbmc.singlet <- FindNeighbors(pbmc.singlet, reduction = "pca", dims = 1:10)
  # pbmc.singlet <- FindClusters(pbmc.singlet, resolution = 0.6, verbose = FALSE)
  # pbmc.singlet <- RunTSNE(pbmc.singlet, reduction = "pca", dims = 1:10)
  # # 
  # if(argv$dimPlot){
  #   plot6 <- DimPlot(pbmc.singlet, group.by = "HTO_classification")
  #   png(paste(graphsPath,"dimPlot.png",sep=""))
  #   print(plot6)
  #   dev.off()
  # }
  
  
}



