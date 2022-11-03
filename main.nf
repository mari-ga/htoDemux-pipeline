#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """\
 Hashtag Demultiplexing - P I P E L I N E
 ===================================
 Input Files:
Seurat:
 RNA-Data: ${params.umi_count}
 HTO-Matrix: ${params.hto_mat}

 Empty Drops - DemuxEM
 HTO-Matrix raw: ${params.hto_raw}

 Hashed Drops:
 HTO-Data: ${params.hashtag_data}
 
 Type of data:
 Raw-data: ${params.rawData}
 demultiplexing: ${params.demultiplexing}
 doublet detection: ${params.doublet_detection}
 filtering data: ${params.cleaning_raw}

 """

include { SEURAT } from './modules/seurat'
include { HASHED_DROPS } from './modules/hashed_drops'
include { DEMUXEM_DEMUL } from './modules/demuxem_demul'
include { HASH_SOLO_DEMUL } from './modules/hash_solo_demul'
include { SOLO_DEMUL } from './modules/solo_demul'
include { EMPTY_DROPS_FLOW } from './modules/empty_drops_flow'
include { ASSIGNMENT_WORKFLOW } from './modules/assignment_flow'
include { CLASSIFICATION_WORKFLOW } from './modules/classification_flow'
include { FINAL_REPORT } from './modules/final_report'

workflow{
  //Params for pre-processing
  rdsObject = Channel.from(params.rdsObject)
  umi = Channel.fromPath(params.umi_count, checkIfExists: true )
  hto_matrix =  Channel.fromPath(params.hto_mat, checkIfExists: true )
  sel_method = Channel.from(params.selection_method)
  ndelim = Channel.from(params.ndelim)
  n_features = Channel.from(params.number_features)
  assay = Channel.from(params.assay)
  a_name = Channel.from(params.assayName)
  margin = Channel.from(params.margin)
  norm_method = Channel.from(params.normalisation_method)
  out_file = Channel.from(params.nameOutputFile)
  objectOutHTO = Channel.from(params.objectOutHTO)
  

  //Params for HTODemux
  quantile_hto = Channel.from(params.quantile_hto)
  kfunc = Channel.from(params.kfunc)
  n_stars = Channel.from(params.nstarts)
  n_samples = Channel.from(params.nsamples)
  seed = Channel.from(params.seed)
  init = Channel.from(params.init)
  out_hto = Channel.from(params.nameOutputFileHTO)
  assignment_hto = Channel.from(params.nameAssignmentFileHTO)

  //Params for MULTI-seq
  quantile_multi = Channel.from(params.quantile_multi)
  autoThresh = Channel.from(params.autoThresh)
  maxIter = Channel.from(params.maxiter)
  qrangeFrom = Channel.from(params.qrangeFrom)
  qrangeTo = Channel.from(params.qrangeTo)
  qrangeBy = Channel.from(params.qrangeBy)
  verbose = Channel.from(params.verbose)
  out_multi = Channel.from(params.nameOutputFileMulti)
  classification_multi = Channel.from(params.nameClassificationFileMulti)

  //Params for HTO-Demul visualisation
  visualisation_seurat = Channel.from(params.visualisationSeurat)
  ridgePlot = Channel.from(params.ridgePlot)
  ridgeNCol = Channel.from(params.ridgeNCol)
  featureScatter = Channel.from(params.featureScatter)
  scatterFeat1 = Channel.from(params.scatterFeat1)
  scatterFeat2 = Channel.from(params.scatterFeat2)
  vlnplot = Channel.from(params.vlnplot)
  vlnFeatures = Channel.from(params.vlnFeatures)
  vlnLog = Channel.from(params.vlnLog)
  tsne = Channel.from(params.tsne)
  tseIdents = Channel.from(params.tseIdents)
  tsneInvert = Channel.from(params.tsneInvert)
  tsneVerbose = Channel.from(params.tsneVerbose)
  tsneApprox = Channel.from(params.tsneApprox)
  tsneDimMax = Channel.from(params.tsneDimMax)
  tsePerplexity = Channel.from(params.tsePerplexity)
  heatmap = Channel.from(params.heatmap)
  heatmapNcells = Channel.from(params.heatmapNcells)

  //Params for Hashed Drops
  hashtag_data = Channel.from(params.hashtag_data)
  nameOutputFileDrops = Channel.from(params.nameOutputFileDrops)
  nameOutputFileHashed = Channel.from(params.nameOutputFileHashed)
  rawData = Channel.from(params.rawData)
  ambient = Channel.from(params.ambient)
  minProp = Channel.from(params.minProp)
  pseudoCount = Channel.from(params.pseudoCount)
  constAmbient = Channel.from(params.constAmbient)
  doubletNmads = Channel.from(params.doubletNmads)
  doubletMin = Channel.from(params.doubletMin)
  confidenMin = Channel.from(params.confidenMin)
  confidentNmads = Channel.from(params.confidentNmads)
  histogram = Channel.from(params.histogram)
  plotLog = Channel.from(params.plotLog)
  empty_drops_result = Channel.from(params.empty_drops_result)
  

  //Params for DemuxEM
  //using filtered HTO matrix - same as Hashed Drops, HtoDemux and Multi-seq
  // using RNA raw matrix - same as empty drops
  alpha = Channel.from(params.alpha)
  alpha_noise = Channel.from(params.alpha_noise)
  tol = Channel.from(params.tol)
  n_threads = Channel.from(params.n_threads)
  min_signal = Channel.from(params.min_signal)
  output_demux = Channel.from(params.output_demux)
  

  //Params for Hash Solo
  hto_data = Channel.from(params.hto_data)
  priors_negative = Channel.from(params.priors_negative)
  priors_singlet = Channel.from(params.priors_singlet)
  priors_doublet = Channel.from(params.priors_doublet)
  output_file = Channel.from(params.output_file)
  output_plot = Channel.from(params.output_plot)

  //params for Solo
  soft = Channel.from(params.soft)
  max_epochs = Channel.from(params.max_epochs)
  lr = Channel.from(params.lr)
  output_solo = Channel.from(params.output_solo)

  //params for Empty Drops
  rna_raw = Channel.from(params.rna_raw)
  hto_raw = Channel.from(params.hto_raw)
  niters = Channel.from(params.niters)
  empty = Channel.from(params.empty)
  lower = Channel.from(params.lower)
  testAmbient = Channel.from(params.testAmbient)
  alpha_empty = Channel.from(params.alpha_empty)
  ignore = Channel.from(params.ignore)
  nameOutputEmpty = Channel.from(params.nameOutputEmpty)
  nameObjectEmpty = Channel.from(params.nameObjectEmpty)
  //general assignment, classification and intermediate file

  output_assignment = Channel.from(params.output_assignment)
  output_classification = Channel.from(params.output_classification)
  output_final =  Channel.from(params.output_final)

  demux = "/Users/mylenemarianagonzalesandre/Development/Results-cluster/Results-Batch1-orig-param/output_demuxEM.csv"


  if(params.demultiplexing == "TRUE"){

    SEURAT(visualisation_seurat,rdsObject,umi,hto_matrix, sel_method,ndelim, n_features, assay, a_name, margin,norm_method,seed, init, out_file, quantile_hto,kfunc, n_stars,n_samples,out_hto,assignment_hto,objectOutHTO,quantile_multi,autoThresh,maxIter,qrangeFrom,qrangeTo,qrangeBy,verbose,out_multi,classification_multi,ridgePlot,ridgeNCol, featureScatter,scatterFeat1,scatterFeat2,vlnplot,vlnFeatures,vlnLog,tsne,tseIdents,tsneInvert,tsneVerbose,tsneApprox,tsneDimMax,tsePerplexity,heatmap,heatmapNcells)
    demux_out_1 = SEURAT.out.HTODEMUX_OUT_1
    demux_out_2 = SEURAT.out.HTODEMUX_OUT_2
    multi_out = SEURAT.out.MULTISEQ_OUT_1

    hashed_drops_out = Channel.empty()
    HASHED_DROPS(empty_drops_result,rawData,hashtag_data,nameOutputFileDrops,nameOutputFileHashed,ambient, minProp,pseudoCount,constAmbient,doubletNmads,doubletMin,confidenMin,confidentNmads,histogram,plotLog,hto_raw,niters,empty,lower,testAmbient,alpha_empty,ignore,nameOutputEmpty,nameObjectEmpty)
    hashed_drops_out = HASHED_DROPS.out.HASHED_DROPS_OUT

    demuxem_out = Channel.empty()
    DEMUXEM_DEMUL(rna_raw,hto_matrix,alpha,alpha_noise,tol,n_threads, min_signal,output_demux)
    demuxem_out = DEMUXEM_DEMUL.out.DEMUXEM_OUT
    

    hash_solo_out = Channel.empty()
    HASH_SOLO_DEMUL(hto_data,priors_negative,priors_singlet,priors_doublet,output_file,output_plot)
    hash_solo_out = HASH_SOLO_DEMUL.out.HASH_SOLO_OUT

    ASSIGNMENT_WORKFLOW(demux_out_2,multi_out,hashed_drops_out,hash_solo_out,demuxem_out,output_assignment)

    CLASSIFICATION_WORKFLOW(demux_out_1,multi_out,hashed_drops_out,hash_solo_out,demuxem_out,output_classification)

  }else{
    print("demultiplexing was not executed")
  }


  if(params.doublet_detection == "TRUE"){
    solo_out = Channel.empty()
    SOLO_DEMUL(umi,soft,max_epochs,lr,output_solo)
    solo_out = SOLO_DEMUL.out.SOLO_OUT
    
  }else{
    print("Solo was not executed")
  }

  if(params.cleaning_raw == "TRUE"){
    EMPTY_DROPS_FLOW(hto_raw,niters,empty,lower,testAmbient,alpha_empty,ignore,nameOutputEmpty,nameObjectEmpty)
  
  }else{
    print("Empty Drops was not executed")
  }

  if (params.demultiplexing == "TRUE" && params.doublet_detection == "TRUE" )
  {
    FINAL_REPORT(CLASSIFICATION_WORKFLOW.out.classification_out,solo_out,output_final)
  }

}