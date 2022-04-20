#!/usr/bin/env nextflow
nextflow.enable.dsl=2
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
log.info """\
 HTODemux - P I P E L I N E
 ===================================
 UMI-Counts: ${params.params}
 HTO-Matrix: ${params.hto_mat}
 """

process htoDemux {
  input:
    path umi_counts
    path hto_matrix
  script:
    """
    Rscript HTODemux-args.R $umi_counts $hto_matrix
    """

}

workflow {
  umi_counts = Channel.fromPath(params.umi_count)
  hto_matrix = Channel.fromPath(params.hto_mat)
}




