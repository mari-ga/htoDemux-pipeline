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
include { ASSIGNMENT_WORKFLOW } from './modules/assignment_flow.nf'
include { EMPTY_DROPS_FLOW } from './modules/empty_drops_flow'

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

  //Params for HTODemul
  quantile_hto = Channel.from(params.quantile_hto)
  kfunc = Channel.from(params.kfunc)
  n_stars = Channel.from(params.nstarts)
  n_samples = Channel.from(params.nsamples)
  seed = Channel.from(params.seed)
  init = Channel.from(params.init)
  out_hto = Channel.from(params.nameOutputFileHTO)

  //Params for MULTI-seq
  quantile_multi = Channel.from(params.quantile_multi)
  autoThresh = Channel.from(params.autoThresh)
  maxIter = Channel.from(params.maxiter)
  qrangeFrom = Channel.from(params.qrangeFrom)
  qrangeTo = Channel.from(params.qrangeTo)
  qrangeBy = Channel.from(params.qrangeBy)
  verbose = Channel.from(params.verbose)
  out_multi = Channel.from(params.nameOutputFileMulti)

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
  rna_data = Channel.from(params.rna_data)
  hto_raw = Channel.from(params.hto_raw)
  alpha = Channel.from(params.alpha)
  alpha_noise = Channel.from(params.alpha_noise)
  tol = Channel.from(params.tol)
  n_threads = Channel.from(params.n_threads)
  min_signal = Channel.from(params.min_signal)
  output_demux = Channel.from(params.output_demux)
  

  //Params for Hash Solo
  hto_data = Channel.from(params.hto_data)
  priors = Channel.from(params.priors)
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
  alpha = Channel.from(params.alpha)
  ignore = Channel.from(params.ignore)
  nameOutputEmpty = Channel.from(params.nameOutputEmpty)
  

if(params.seurat == 'TRUE'){
  seurat = SEURAT(umi,hto_matrix, sel_method,ndelim, n_features, assay, a_name, margin,norm_method,seed, init, out_file, quantile_hto,kfunc, n_stars,n_samples,out_hto,quantile_multi,autoThresh,maxIter,qrangeFrom,qrangeTo,qrangeBy,verbose,out_multi, ridgePlot,ridgeNCol, featureScatter,scatterFeat1,scatterFeat2,vlnplot,vlnFeatures,vlnLog,tsne,tseIdents,tsneInvert,tsneVerbose,tsneApprox,tsneDimMax,tsePerplexity,heatmap,heatmapNcells)
}

if(params.hashedMode == 'TRUE'){
  hashed = HASHED_DROPS(umi,hto_matrix,nameOutputFileDrops,nameOutputFileHashed,ambient, minProp,pseudoCount,constAmbient,doubletNmads,doubletMin,confidenMin,confidentNmads,combinations,histogram,plotLog,empty,lower,testAmbient,nameOutputEmpty)
}

if(params.demuxem_mode == 'TRUE'){
  demux = DEMUXEM_DEMUL(rna_data,hto_raw,alpha,alpha_noise,tol,n_threads, min_signal,output_demux)
}

if(params.hash_solo_mode == 'TRUE'){
  hash_solo = HASH_SOLO_DEMUL(hto_data,priors,output_file,output_plot)
}

if(params.solo_mode == 'TRUE'){
  solo SOLO_DEMUL(umi,soft,max_epochs,lr,output_solo)
}

if(params.empty_drops_mode == "TRUE")
{
  EMPTY_DROPS_FLOW(rna_raw,niters,empty,lower,testAmbient,alpha,ignore,nameOutputEmpty)
}


if(params.general_assignment == 'TRUE'){
  ASSIGNMENT_WORKFLOW(SEURAT.out[0],SEURAT.out[4],HASHED_DROPS.out[1],DEMUXEM_DEMUL.out,HASH_SOLO_DEMUL.out[0],SOLO_DEMUL.out)

}


}


//params.outdir = '/home/icb/mariana.gonzales/pipeline/demultiplex-pipeline/results/'