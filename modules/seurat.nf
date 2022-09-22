params.outdir = 'results'
//both, hto, multi


include { PREPROCESS } from './preprocess'
include { HTODEMUL } from './htodemul'
include { HTO_VISUALISATION } from './hto_visualisation'
include { MULTI_SEQ } from './multi_complete'


workflow SEURAT{
    take:
        umi_matrix
        hto_matrix
        sel_method
        ndelim
        n_features
        assay 
        a_name 
        margin 
        norm_method 
        seed
        init
        out_file

        quantile_hto
        kfunc
        n_stars 
        n_samples
        out_hto

        quantile_multi
        autoThresh
        maxIter
        qrangeFrom
        qrangeTo
        qrangeBy
        verbose
        out_multi

        ridgePlot
        ridgeNCol
        featureScatter 
        scatterFeat1 
        scatterFeat2 
        vlnplot
        vlnFeatures 
        vlnLog 
        tsne 
        tseIdents 
        tsneInvert 
        tsneVerbose 
        tsneApprox 
        tsneDimMax 
        tsePerplexity 
        heatmap 
        heatmapNcells 


    main:
    PREPROCESS(umi_matrix, hto_matrix,ndelim,sel_method, n_features,assay,a_name, margin, norm_method,out_file)

  
      
    HTODEMUL(PREPROCESS.out,quantile_hto,kfunc,n_stars,n_samples,seed, init,out_hto)

    
    if(params.visualisationSeurat == 'TRUE')
      {
        HTO_VISUALISATION(HTODEMUL.out[0],a_name, ridgePlot,ridgeNCol,featureScatter,scatterFeat1,scatterFeat2,vlnplot,vlnFeatures,vlnLog,tsne,tseIdents,tsneInvert,tsneVerbose,tsneApprox,tsneDimMax,tsePerplexity,heatmap,heatmapNcells)
      }

    MULTI_SEQ(umi_matrix, hto_matrix,ndelim,sel_method, n_features,assay,a_name, margin, norm_method,quantile_multi, autoThresh,  maxIter,qrangeFrom,qrangeTo,qrangeBy,verbose,out_multi)
    
        
  

    
}