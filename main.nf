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
 Intermediate object name: ${params.nameOutputFile}
 Intermediate object path: ${params.demulOutPath}
 Results: ${params.outdir}
 HTOOut: ${params.htoDemuxOutPath}
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
    path 'object.rds'  
  
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
  Rscript $baseDir/R/pre_processing.R ${umiFile} ${htoFile} ${selectMethod} ${numberFeatures} ${assay} ${assayName} ${margin} ${normalisationMethod} ${demulOutPath} ${fileName} > object.rds
  
  """
}


process htoDemux{
  publishDir path: "$params.outdir"
  input:
    path preprocess_object
    val quantile_hto
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
    def objectFile_hto = "--seuratObjectPath $preprocess_object"
    def quantile_hto = "--quantile $quantile_hto"
    
    def kfunc = "--kfunc  $kfunc"
    def nstarts = "--nstarts $nstarts"
    def nsamples = "--nsamples $nsamples " 
    def htoOutpath = "--htoDemuxOutPath $htoDemuxOutPath"
    def nameFileHTO = "--nameOutputFileHTO $nameOutputFileHTO "
    


    """
      Rscript $baseDir/R/HTODemux-args.R ${preprocess_object} ${quantile_hto} ${kfunc} ${nstarts} ${nsamples} ${htoOutpath} ${nameFileHTO}
    """
}

process multiSeq{
  publishDir path: "$params.outdir"
  input:
  path preprocess_object
  val quantile_multi
  //For Multi-seq
  val autoThresh
  val maxiter
  val qrangeFrom
  val qrangeTo
  val qrangeBy
  val verbose
  val multiSeqOutPath
  val nameOutputFileMulti

  def objectFile_multi = "--seuratObjectPath $preprocess_object"
  def quantile_multi = "--quantile $quantile"
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
    Rscript $baseDir/R/MULTI-seq.R ${objectFile_multi} ${quantile_multi} ${autoThresh} ${maxiter} ${qrangeFrom} ${qrangeTo} ${qrangeBy} ${verbose} ${multiSeqOutPath} ${nameOutputFileMulti}
  """

}


process show{
  input:
  path preprocess_object


  script:
  """
  echo recibido: $preprocess_object
  
  """


}

//Subworkflows
workflow pre_processing{
    def umi = Channel.fromPath(params.umi_count, checkIfExists: true )
    def hto_matrix =  Channel.fromPath(params.hto_mat, checkIfExists: true )
  main:
    preProcess(umi,hto_matrix,params.selection_method, params.number_features, params.assay, params.assayName, params.margin, params.normalisation_method, params.demulOutPath, params.nameOutputFile)
  emit:
	  output_object = preProcess.out.view({ "Created: $it" })
}


// workflow demul_htoDemux{
//  take: 
//     path pre_processing.out.output_object
//   main:
//       htoDemux(pre_processing.out.output_object,params.quantile_hto, params.kfunc, params.nstarts, params.nsamples, params.htoOutpath, params.nameFileHTO)
//   emit:
//     htoDemux_out_rds = htoDemux.out[0]
//     htoDemux_out_csv = htoDemux.out[1]

// }

// workflow demul_multiSeq{
//   main:
//       multiSeq(pre_processing.out.output_object, params.quantile_multi, params.autoThresh, params.maxiter, params.qrangeFrom, params.qrangeTo, params.qrangeBy, params.verbose, params.multiSeqOutPath, params.nameOutputFileMulti)
//   emit:
//     multiSeq_out_rds = multiSeq.out[0]
//     multiSeq_out_csv = multiSeq.out[1]
// }



//Main workflow
workflow{
  // take:

  main:
  out_pre = pre_processing()
  //pre_processing.out.view({ "Received: $it" })
  myFileChannel = Channel.fromPath(out_pre)
  println(myFileChannel)
  //htoDemux(myFileChannel,params.quantile_hto, params.kfunc, params.nstarts, params.nsamples, params.htoDemuxOutPath, params.nameOutputFileHTO)
 
  //   demul_htoDemux(out_pre)
    //if mode on -> cual (ifs anidados)
   // if ()

   // else()




}


workflow.onComplete { 
  println ("Done")
}



