#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
 * Input: 
    --umi-matrix.rds <path>
    --hashtag-counts-matrix.rds <path>
    Parameters for HTODemux process:
    --selection_method
    --number_features
    --normalisation_method
    --margin
    --assay
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
    val assayName
    var kfunc
    path htoDemuxOutPath
    path graphs
    path nameOutputFile
    val nameOutputFile
    val ridgePlot
    val ridgeNCol
    val featureScatter
    val scatterFeat1
    val scatterFeat2

  output:
    stdout seurat_object
  script:
    """
    Rscript HTODemux-args.R --fileUmi ${umi_counts} --fileHto ${hto_matrix}  --selectMethod ${selection_method} --numberFeatures ${number_features} --normalisationMethod ${normalisation_method} --margin ${margin} 
    --assay ${assay} --assayName ${assayName} --kfunc ${kfunc} --htoDemuxOutPath ${htoDemuxOutPath} --graphs ${graphs} --nameOutputFile ${nameOutputFile} --ridgePlot ${ridgePlot} --ridgeNCol ${ridgeNCol}
    --featureScatter ${featureScatter} --scatterFeat1 ${scatterFeat1} --scatterFeat2 ${scatterFeat2}
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
  assayName = channel.value(params.assayName)
  kfunc = channel.value(params.kfunc)
  htoDemuxOutPath = Channel.fromPath(params.htoDemuxOutPath)
  graphs = Channel.fromPath(params.graphs)
  nameOutputFile = channel(params.nameOutputFile)
  ridgePlot = channel(params.ridgePlot)
  ridgeNCol = channel(params.ridgeNCol)
  featureScatter = channel(params.featureScatter)
  scatterFeat1 = channel(params.scatterFeat1)
  scatterFeat2 = channel(params.scatterFeat2)
}




