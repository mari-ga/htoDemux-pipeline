#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
 * Input: 
    --umi-matrix.rds <path>
    --hashtag-counts-matrix.rds <path>
    Parameters for Pre-processing (Seurat object creation):
    --selection_method
    --number_features
    --normalisation_method
    --margin
    --assay
    --assayName
    --demulOutPath
    --nameOutputFile
    
    
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
    val nameOutputFile

  output:
    path 'object' 
  
  script:
  def umiFile = "--fileUmi $umi_counts"
  def htoFile = "--fileHto $hto_matrix"
  def selectMethod = "--selectMethod $selection_method"
  def numberFeatures = "--numberFeatures $number_features"
  def assay = "--assay $assay"
  def assayName = "--assayName $assayName"
  def margin = "--margin $margin"
  def normalisationMethod = "--normalisationMethod $normalisation_method"
  def demulOutPath = "--demulOutPath $demulOutPath"
  def fileName = " --nameOutputFile $nameOutputFile"

  """
Rscript $baseDir/pre_processing.R ${umiFile} ${htoFile} ${selectMethod} ${numberFeatures} ${assay} ${assayName} ${margin} ${normalisationMethod} ${demulOutPath} ${fileName}
  """

}




workflow{
    def umi = Channel.fromPath(params.umi_count)
    def hto_matrix =  Channel.fromPath(params.hto_mat)
  main:
    preProcess(umi,hto_matrix,params.selection_method, params.number_features, params.assay, params.assayName, params.margin, params.normalisation_method, params.demulOutPath, params.nameOutputFile)
  emit:
	    preProcess.out
}


workflow.onComplete { 
  println ("Done")
}



