#!/usr/bin/env nextflow
nextflow.enable.dsl=2

num = Channel.from( 1, 2, 3 )

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

