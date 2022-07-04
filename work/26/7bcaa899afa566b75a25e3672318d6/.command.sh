#!/bin/bash -ue
Rscript /Users/mylenemarianagonzalesandre/Development/Bachelor-Thesis/demultiplex-pipeline/R/pre_processing.R --fileUmi pbmc_umi_mtx.rds --fileHto pbmc_hto_mtx.rds --selectMethod mean.var.plot --numberFeatures 2000 --assay HTO --assayName HTO --margin 2 --normalisationMethod CLR --nameOutputFile object
