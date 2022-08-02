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

workflow{
  //Params for pre-processing
  umi = Channel.fromPath(params.umi_count, checkIfExists: true )
  hto_matrix =  Channel.fromPath(params.hto_mat, checkIfExists: true )
  sel_method = Channel.from(params.selection_method)
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

  //Params for DemuxEM
  rna_data = Channel.from(params.rna_data)
  hto_mat_em = Channel.from(params.hto_mat_em)
  threads = Channel.from(params.threads)
  alpha = Channel.from(params.alpha)
  alpha_noise = Channel.from(params.alpha_noise)
  min_signal = Channel.from(params.min_signal)
  tol = Channel.from(params.tol)


  //incluir un if - seurat entra a todo el mambo
  SEURAT(umi,hto_matrix, sel_method, n_features, assay, a_name, margin,norm_method, out_file, quantile_hto,kfunc, n_stars,n_samples,out_hto,quantile_multi,autoThresh,maxIter,qrangeFrom,qrangeTo,qrangeBy,verbose,out_multi, ridgePlot,ridgeNCol, featureScatter,scatterFeat1,scatterFeat2,vlnplot,vlnFeatures,vlnLog,tsne,tseIdents,tsneInvert,tsneVerbose,tsneApprox,tsneDimMax,tsePerplexity,heatmap,heatmapNcells)

  HASHED_DROPS(umi,hto_matrix,nameOutputFileDrops,nameOutputFileHashed,ambient, minProp,pseudoCount,constAmbient,doubletNmads,doubletMin,confidenMin,confidentNmads,combinations)

  //DEMUXEM(rna_data,hto_mat_em,min_signal,alpha,alpha_noise,tol,threads,)
}
