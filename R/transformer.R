#!/usr/bin/env Rscript

library(argparser, quietly=TRUE)


#This script is a transformer for the UMI and HTO counts matrices in .rds to csv and hdf5 respectively
p <- arg_parser("Parameters transformer")