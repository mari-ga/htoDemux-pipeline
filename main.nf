#!/usr/bin/env nextflow
nextflow.enable.dsl=2


"""
 * Input: 
    --umi-matrix.rds
    --hashtag-counts-matrix.rds

  * Output    
    -- 
"""
umi_chanel = Channel.fromPath(params.umi_count)
hto_chanel = Channel.fromPath(params.htos_mat)

input:

file umi_counts from umi_chanel

output:
file 'table_classification' 

process hto-demultiplex{

echo "Working"

}
