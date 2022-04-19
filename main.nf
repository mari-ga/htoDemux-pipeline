#!/usr/bin/env nextflow
nextflow.enable.dsl=2


"""
 * Input: 
    --umi-matrix.rds <path>
    --hashtag-counts-matrix.rds <path>
    --help

  * Output    
    --outdir_binary
    --outdir_ascii
    --out_stdout
"""
umi_chanel = Channel.fromPath(params.umi_count)
hto_chanel = Channel.fromPath(params.htos_mat)


process hto-demultiplex{

input:

    path umi_counts from umi_chanel
    path hto_matrix from hto_chanel

output:

    //Binary
    file 'table_classification' 
    """
        Allow different types of outputs:
            * Binary file is possible from S4 class
            * Ascii file also possible from S4
    """
script:
    """
    connect to R script as file?
    """
echo "Working"

}

process hto-visualisation{


    
}
