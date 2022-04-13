#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.umi = "/Users/mylenemarianagonzalesandre/Development/Bachelor-Thesis/R/Week-1/HTODemuxFiles/pbmc_umi_mtx.rds"

"""
 * Input: 
    * UMI-matrix.rds
    * hashtag-counts-matrix.rds
"""
umi_chanel = Channel.fromPath(  )

input:

file umi_counts from 


process hto-setup{


}

process hto-independent-assay{


}

process hto-demultiplex{


}
process hto-visualise{


}