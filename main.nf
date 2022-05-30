#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
 * Input: 
    Files for Pre Processing:
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
    Files for Seurat Demultiplexing:
    --seuratObjectPath (output from pre processing)



    
    
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
    path 'object', emit: preprocess_object
  
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


process htoDemux{
  input:
    path object
    val quantile
    //For HTODemux
    val kfunc
    val nstarts
    val nsamples
    val htoDemuxOutPath
    val nameOutputFileHTO

  output:
    path 'rds_result_hto'
    path 'csv_result_hto'
  
  script:
  //Parameters in common amongst both algoriths
    def objectFile = "--seuratObjectPath $object"
    def quantile = "--quantile $quantile"
    
    def kfunc = "--kfunc  $kfunc"
    def nstarts = "--nstarts $nstarts"
    def nsamples = "--nsamples $nsamples " 
    def htoOutpath = "--htoDemuxOutPath $htoDemuxOutPath"
    def nameFileHTO = "--nameOutputFileHTO $nameOutputFileHTO "
    
    """
      echo 'Running HTODemux'
      Rscript $baseDir/R/HTODemux-args.R ${objectFile} ${quantile} ${nstarts} ${nsamples} ${htoOutpath} ${nameFileHTO}
    """
}

process multiSeq{
  path object
  val quantile
  //For Multi-seq
  val autoThresh
  val maxiter
  val qrangeFrom
  val qrangeTo
  val qrangeBy
  val verbose
  val multiSeqOutPath
  val nameOutputFileMulti

  def objectFile = "--seuratObjectPath $object"
  def quantile = "--quantile $quantile"
  def objectFile = "--seuratObjectPath $object"
  def quantile = "--quantile $quantile"
  def autoThresh = "--autoThresh $autoThresh"
  def maxiter = "--maxiter $maxiter"
  def qrangeFrom = "--qrangeFrom $qrangeFrom"
  def qrangeTo = "--qrangeTo $qrangeTo"
  def qrangeBy = "--qrangeBy $qrangeBy"
  def verbose = "--verbose $verbose"
  def multiSeqOutPath = "--multiSeqOutPath $multiSeqOutPath"
  def nameOutputFileMulti = "--nameOutputFileMulti $nameOutputFileMulti"


  output:
    path 'rds_result_multi'
    path 'csv_result_multi'

  script:
  """
    echo 'Running MULTI-seq'
    Rscript $baseDir/R/MULTI-seq.R ${objectFile} ${quantile} ${autoThresh} ${maxiter} ${qrangeFrom} ${qrangeTo} ${qrangeBy} ${verbose} ${multiSeqOutPath} ${nameOutputFileMulti}
  """

}



//Subworkflows
workflow pre_processing{
    def umi = Channel.fromPath(params.umi_count)
    def hto_matrix =  Channel.fromPath(params.hto_mat)
  main:
    preProcess(umi,hto_matrix,params.selection_method, params.number_features, params.assay, params.assayName, params.margin, params.normalisation_method, params.demulOutPath, params.nameOutputFile)
  emit:
	    output_object = preProcess.out
}


workflow demul_htoDemux{
  main:
      htoDemux(pre_processing.out.output_object,params.quantile, params.kfunc, params.nstarts, params.nsamples, params.htoOutpath, params.nameFileHTO)
  emit:
    htoDemux_out_rds = htoDemux.out[0]
    htoDemux_out_csv = htoDemux.out[1]

}

workflow demul_multiSeq{
  main:
      multiSeq(pre_processing.out.output_object, params.quantile, params.autoThresh, params.maxiter, params.qrangeFrom, params.qrangeTo, params.qrangeBy, params.verbose, params.multiSeqOutPath, params.nameOutputFileMulti)
  emit:
    multiSeq_out_rds = multiSeq.out[0]
    multiSeq_out_csv = multiSeq.out[1]
}



//Main workflow
workflow{

  main:
    pre_processing()
    //if mode on -> cual (ifs anidados)
    if ()

    else()

}


workflow.onComplete { 
  println ("Done")
}



