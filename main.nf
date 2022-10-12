#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """\
 Hashtag Demultiplexing - P I P E L I N E
 ===================================
 Input Files:

 UMI-Counts: ${params.umi_count}
 HTO-Matrix: ${params.hto_mat}

 """

include { SEURAT } from './modules/seurat'
include { HASHED_DROPS } from './modules/hashed_drops'
include { DEMUXEM_DEMUL } from './modules/demuxem_demul'
include { HASH_SOLO_DEMUL } from './modules/hash_solo_demul'
include { SOLO_DEMUL } from './modules/solo_demul'
include { ASSIGNMENT_WORKFLOW } from './modules/assignment_flow'
include { EMPTY_DROPS_FLOW } from './modules/empty_drops_flow'
include { EMPTY_REPORT } from './modules/report_maker'
include { CLASSIFICATION_WORKFLOW } from './modules/classification_flow'

workflow{
  //Params for pre-processing
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
  nameOutputFileDrops = Channel.from(params.nameOutputFileDrops)
  nameOutputFileHashed = Channel.from(params.nameOutputFileHashed)
  ambient = Channel.from(params.ambient)
  minProp = Channel.from(params.minProp)
  pseudoCount = Channel.from(params.pseudoCount)
  constAmbient = Channel.from(params.constAmbient)
  doubletNmads = Channel.from(params.doubletNmads)
  doubletMin = Channel.from(params.doubletMin)
  confidenMin = Channel.from(params.confidenMin)
  confidentNmads = Channel.from(params.confidentNmads)
  combinations = Channel.from(params.combinations)
  histogram = Channel.from(params.histogram)
  plotLog = Channel.from(params.plotLog)
  

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
  empty_drops_mode = Channel.from(params.empty_drops_mode)
  niters = Channel.from(params.niters)
  empty = Channel.from(params.empty)
  lower = Channel.from(params.lower)
  testAmbient = Channel.from(params.testAmbient)
  alpha_empty = Channel.from(params.alpha_empty)
  ignore = Channel.from(params.ignore)
  nameOutputEmpty = Channel.from(params.nameOutputEmpty)

  empty_col = Channel.from(params.col_1)
  output_assignment = Channel.from(params.output_assignment)
  output_classification = Channel.from(params.output_classification)
//The next lines correspond to the workflows for all the tools in the project
//It is not obligatory to use all of them at once

if(params.seurat == 'TRUE'){
  SEURAT(umi,hto_matrix, sel_method,ndelim, n_features, assay, a_name, margin,norm_method,seed, init, out_file, quantile_hto,kfunc, n_stars,n_samples,out_hto,assignment_hto,objectOutHTO,quantile_multi,autoThresh,maxIter,qrangeFrom,qrangeTo,qrangeBy,verbose,out_multi,classification_multi,ridgePlot,ridgeNCol, featureScatter,scatterFeat1,scatterFeat2,vlnplot,vlnFeatures,vlnLog,tsne,tseIdents,tsneInvert,tsneVerbose,tsneApprox,tsneDimMax,tsePerplexity,heatmap,heatmapNcells)
}

if(params.hashedMode == 'TRUE'){
  HASHED_DROPS(hto_matrix,nameOutputFileDrops,nameOutputFileHashed,ambient, minProp,pseudoCount,constAmbient,doubletNmads,doubletMin,confidenMin,confidentNmads,combinations,histogram,plotLog)
}

if(params.demuxem_mode == 'TRUE'){
  DEMUXEM_DEMUL(rna_raw,hto_matrix,alpha,alpha_noise,tol,n_threads, min_signal,output_demux)
}

if(params.hash_solo_mode == 'TRUE'){
  HASH_SOLO_DEMUL(hto_data,priors_negative,priors_singlet,priors_doublet,output_file,output_plot)
}

if(params.solo_mode == 'TRUE'){
  SOLO_DEMUL(umi,soft,max_epochs,lr,output_solo)
}

if(params.empty_drops_mode == "TRUE")
{
  EMPTY_DROPS_FLOW(rna_raw,niters,empty,lower,testAmbient,alpha_empty,ignore,nameOutputEmpty)
}




ASSIGNMENT_WORKFLOW(SEURAT.out.HTODEMUX_OUT_2, SEURAT.out.MULTISEQ_OUT_1,HASHED_DROPS.out.HASHED_DROPS_OUT,HASH_SOLO_DEMUL.out.HASH_SOLO_OUT,output_assignment)

CLASSIFICATION_WORKFLOW(SEURAT.out.HTODEMUX_OUT_1,SEURAT.out.MULTISEQ_OUT_1,HASHED_DROPS.out.HASHED_DROPS_OUT,HASH_SOLO_DEMUL.out.HASH_SOLO_OUT,output_classification)

}


//params.outdir = '/home/icb/mariana.gonzales/pipeline/demultiplex-pipeline/results/'