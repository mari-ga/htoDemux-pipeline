#!/usr/bin/env nextflow
nextflow.enable.dsl=2


"""
 * Input: 
    * UMI-matrix.rds
    * hashtag-counts-matrix.rds
"""
umi_chanel = Channel.fromPath(params.umi_counts)

input:

file umi_counts from umi_chanel

output:
file 'table_classification' 

process hto-demultiplex{

echo "Working"

}
