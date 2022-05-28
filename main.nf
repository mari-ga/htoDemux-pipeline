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
 Hashtag Demultiplexing - P I P E L I N E
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
    path 'object' into seurat_demul
  
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
  Rscript $baseDir/R/pre_processing.R ${umiFile} ${htoFile} ${selectMethod} ${numberFeatures} ${assay} ${assayName} ${margin} ${normalisationMethod} ${demulOutPath} ${fileName}
  """
}


process demultiplexing
{
  input:
    path 'object' from seurat_demul
    val mode 
    val quantile
    //For HTODemux
    val kfunc
    val nstarts
    val nsamples
    val htoDemuxOutPath
    val nameOutputFileHTO

    //For Multi-seq
    val autoThresh
    val maxiter
    val qrangeFrom
    val qrangeTo
    val qrangeBy
    val verbose
    val multiSeqOutPath
    val nameOutputFileMulti

  output:
    path 'rds_result'
    path 'csv_result'
  
  script:
  //Parameters in common amongst both algoriths
    def objectFile = "--seuratObjectPath $object"
    def quantile = "--quantile $quantile"
    
    def kfunc = "--kfunc  $kfunc"
    def nstarts = "--nstarts $nstarts"
    def nsamples = "--nsamples $nsamples " 
    def htoOutpath = "--htoDemuxOutPath $htoDemuxOutPath"
    def nameFileHTO = "--nameOutputFileHTO $nameOutputFileHTO "
    
    def autoThresh = "--autoThresh $autoThresh"
    def maxiter = "--maxiter $maxiter"
    def qrangeFrom = "--qrangeFrom $qrangeFrom"
    def qrangeTo = "--qrangeTo $qrangeTo"
    def qrangeBy = "--qrangeBy $qrangeBy"
    def verbose = "--verbose $verbose"
    def multiSeqOutPath = "--multiSeqOutPath $multiSeqOutPath"
    def nameOutputFileMulti = "--nameOutputFileMulti $nameOutputFileMulti"


    if(mode == "htodemux")
      """
        echo 'Running HTODemux'
        Rscript $baseDir/R/HTODemux-args.R ${objectFile} ${quantile} ${nstarts} ${nsamples} ${htoOutpath} ${nameFileHTO}
      """

    else if(mode == "multiseq")
      """
        echo 'Running MULTI-seq'
        Rscript $baseDir/R/MULTI-seq.R ${objectFile} ${quantile} ${autoThresh} ${maxiter} ${qrangeFrom} ${qrangeTo} ${qrangeBy} ${verbose} ${multiSeqOutPath} ${nameOutputFileMulti}
      """

    else if(mode == "both")
      """
        echo 'Running HTODemux'
        Rscript $baseDir/R/HTODemux-args.R ${objectFile} ${quantile} ${nstarts} ${nsamples} ${htoOutpath} ${nameFileHTO}
        echo 'Running MULTI-seq'
        Rscript $baseDir/R/MULTI-seq.R ${objectFile} ${quantile} ${autoThresh} ${maxiter} ${qrangeFrom} ${qrangeTo} ${qrangeBy} ${verbose} ${multiSeqOutPath} ${nameOutputFileMulti}
      """

    else
          error "Invalid alignment mode: ${mode}"
  


}


//Subworkflows
workflow pre_processing{
    def umi = Channel.fromPath(params.umi_count)
    def hto_matrix =  Channel.fromPath(params.hto_mat)
  main:
    preProcess(umi,hto_matrix,params.selection_method, params.number_features, params.assay, params.assayName, params.margin, params.normalisation_method, params.demulOutPath, params.nameOutputFile)
  emit:
	    preProcess.out
}
workflow demul_seurat{


}



//Main workflow
workflow{
   take:
    def umi = Channel.fromPath(params.umi_count)
    def hto_matrix =  Channel.fromPath(params.hto_mat)
  main:
    pre_processing(umi,hto_matrix)

}


workflow.onComplete { 
  println ("Done")
}



