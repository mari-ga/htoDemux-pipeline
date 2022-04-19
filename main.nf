#!/usr/bin/env nextflow
nextflow.enable.dsl=2

umi_chanel = Channel.fromPath(params.umi_count)
hto_chanel = Channel.fromPath(params.htos_mat)

/*
 * Input: 
    --umi-matrix.rds <path>
    --hashtag-counts-matrix.rds <path>
    --help

  * Output    
    --outdir_binary
    --outdir_ascii
    --out_stdout
*/

process htoDemultiplex{
 input:
  val x from num

  "echo process job $x"


}




input:

    path umi_counts from umi_chanel
    path hto_matrix from hto_chanel

output:

    //Binary
    file 'table_classification' 
    /*
        Allow different types of outputs:
            * Binary file is possible from S4 class
            * Ascii file also possible from S4
    */

process hto-visualisation{



}


