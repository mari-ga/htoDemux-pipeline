#!/bin/bash -ue
Rscript /Users/mylenemarianagonzalesandre/Development/Bachelor-Thesis/demultiplex-pipeline/R/HTODemux-args.R --seuratObjectPath object.rds --quantile 0.99 --kfunc  kmeans --nstarts 100 --nsamples 100 --nameOutputFileHTO resultHTO
