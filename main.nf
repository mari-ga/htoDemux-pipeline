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
    val selection_method
    val number_features
    val normalisation_method
    val margin
    val assay

  output:
    stdout seurat_object
  script:
    """
    Rscript HTODemux-args.R ${umi_counts} ${hto_matrix} ${selection_method} ${number_features} ${normalisation_method} ${margin} ${assay}
    """

}

workflow {
  umi_counts = Channel.fromPath(params.umi_count)
  hto_matrix = Channel.fromPath(params.hto_mat)
  selection_method = channel.value(params.selection_method)
  number_features = channel.value(params.number_features)
  normalisation_method = channel.value(params.normalisation_method)
  margin = channel.value(params.margin)
  assay = channel.value(params.assay)
}




