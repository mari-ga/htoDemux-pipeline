
#install.packages("KernSmooth",repos=("http://cran.rstudio.com"))
#install.packages("reshape2",repos=("http://cran.rstudio.com"))
#install.packages("rTSNE",repos=("http://cran.rstudio.com"))
#install.packages("stringdist",repos=("http://cran.rstudio.com"))
#install.packages("ShortRead",repos=("http://cran.rstudio.com"))
devtools::install_github('chris-mcginnis-ucsf/MULTI-seq')
install.packages('devtools')
library(deMULTIplex)

#Inputs
#Bar file
#Cell file

#Outputs
#Graph -> TSNE
#Graph -> 2nd TSNE
#Graph -> results 




##---------------- Section 1 - Input files -----------------
## Define vectors for reference barcode sequences and cell IDs

bar.ref <- load("/Users/mylenemarianagonzalesandre/Development/Bachelor-Thesis/nextflow-files/htoDemux-pipeline/data/multi-seq/LMOlist.csv")
cell.id.vec <- load("/Users/mylenemarianagonzalesandre/Development/Bachelor-Thesis/nextflow-files/htoDemux-pipeline/data/multi-seq/cellIDs.txt")



readTable <- MULTIseq.preProcess(R1 = '/Users/mylenemarianagonzalesandre/Development/Bachelor-Thesis/nextflow-files/htoDemux-pipeline/data/multi-seq/ACAGTG_S3_L001_R1_001.fastq.gz', R2 = '/Users/mylenemarianagonzalesandre/Development/Bachelor-Thesis/nextflow-files/htoDemux-pipeline/data/multi-seq/ACAGTG_S3_L001_R2_001.fastq.gz', cellIDs = cell.id.vec, cell=c(1,16), umi=c(17,28), tag=c(1,8))

##---------------- Section 2 - Inspect Barcode quality -----------------

## Visualize barcode space
bar.tsne <- barTSNE(bar.table[,1:96]) 
## Note: Exclude columns 97:98 (assuming 96 barcodes were used) which provide total barcode UMI counts for each cell. 

pdf("bc.check.pdf")
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}
dev.off()


##---------------- Section 3 - Sample classification -----------------

## Round 1 -----------------------------------------------------------------------------------------------------
## Perform Quantile Sweep
bar.table.full <- bar.table[,1:96]
good.bars <- paste("Bar",1:90,sep="")  # NOTE: In this hypothetical example, barcodes 91-96 were not detected
bar.table <- bar.table.full[, good.bars]  # Remove missing bars and summary columns
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))


## Finalize round 1 classifications, remove negative cells
round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

## Round 2 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
round2.calls <- classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema))
neg.cells <- c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])

## Repeat until all no negative cells remain (usually 3 rounds)...
final.calls <- c(round2.calls, rep("Negative",length(neg.cells)))
names(final.calls) <- c(names(round2.calls),neg.cells)


##---------------- Section 4 - Semi-supervised negative cell reclassification  -----------------

## Perform semi-supervised negative cell reclassification
reclass.cells <- findReclassCells(bar.table.full, names(final.calls)[which(final.calls=="Negative")])
reclass.res <- rescueCells(bar.table.full, final.calls, reclass.cells)


## Visualize Results
ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) + 
  geom_point() + xlim(c(nrow(reclass.res)-1,1)) + 
  ylim(c(0,1.05)) +
  geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2)

## Finalize negative cell rescue results
final.calls.rescued <- final.calls
rescue.ind <- which(reclass.cells$ClassStability >= 16) ## Note: Value will be dataset-specific
final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]
