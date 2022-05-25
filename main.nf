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
    --val assayName
    --kfunc
    --htoDemuxOutPath
    --graphs
    --nameOutputFile
    --nameOutputFile
    --ridgePlot
    --ridgeNCol
    --featureScatter
    --scatterFeat1
    --scatterFeat2
    --vlnplot
    --vlnFeatures
    --vlnLog
    --tsne
    --tseIdents
    --tsneInvert
    --tsePerplexity
    --heatmap
    --heatmapNcells
    --cluster
    --clusterIdents
    --clusterSelMethod
    --reductionMethod
    --reductionDims
    --reductionResol
    --dimPlot

  * Output    
    --outdir_binary
    --outdir_ascii
    --out_stdout
*/
log.info """\
 Demultiplexing - P I P E L I N E
 ===================================
 Input Files:

 UMI-Counts: ${params.umi_count}
 HTO-Matrix: ${params.hto_mat}
 """

process preProcess{
  input:
    path umi_counts
    path hto_matrix
    val selection_method
    val number_features
    val assay
    val assayName
    val margin
    val normalisation_method
    val demulOutPath
    var nameOutputFile

  output:
    path 'object' 
  
  script:

  """
    Rscript pre-processing.R --fileUmi ${umi_counts} --fileHto ${hto_matrix} --selectMethod ${selection_method}
     --numberFeatures ${number_features} --assay ${assay} --assayName ${assayName}  --margin ${margin} 
     --normalisationMethod ${normalisation_method} --demulOutPath ${demulOutPath} --nameOutputFile ${nameOutputFile}
  """

}




workflow {
  umi = Channel.fromPath(params.umi_count)
  hto_matrix =  Channel.fromPath(params.hto_mat)
  selection_method = params.selection_method
  number_features = params.number_features
  normalisation_method = params.normalisation_method
  margin = params.margin
  assay = params.assay
  assayName = params.assayName
  demulOutPath = params.demulOutPath
  nameOutputFile = params.nameOutputFile

  preProcess(umi,hto_matrix, selection_method, number_features, assay, assayName, margin, normalisation_method,demulOutPath, nameOutputFile)
  
}


workflow.onComplete { 
  println ("Done")
}



