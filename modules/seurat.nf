params.outdir = 'results'
//both, hto, multi


include { PREPROCESS } from './preprocess'
include { HTODEMUL } from './htodemul'
include { MULTISEQ } from './multiseq'
include { HTO_VISUALISATION } from './hto_visualisation'

workflow SEURAT{
    take:
        umi_matrix
        hto_matrix
        sel_method
        n_features
        assay 
        a_name 
        margin 
        norm_method 
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
    PREPROCESS(umi_matrix, hto_matrix,sel_method, n_features,assay,a_name, margin, norm_method,out_file)
    if( params.mode == 'multi' )
    {
       MULTISEQ(PREPROCESS.out, quantile_multi, autoThresh,  maxIter,qrangeFrom,qrangeTo,qrangeBy,verbose,out_multi )
   
    }
        
    else {
        HTODEMUL(PREPROCESS.out,quantile_hto,kfunc,n_stars,n_samples,out_hto)
        if(params.visualisationSeurat == 'TRUE')
        {
            HTO_VISUALISATION(HTODEMUL.out[0],ridgePlot,ridgeNCol,featureScatter,scatterFeat1,scatterFeat2,vlnplot,vlnFeatures,vlnLog,tsne,tseIdents,tsneInvert,tsneVerbose,tsneApprox,tsneDimMax,tsePerplexity,heatmap,heatmapNcells)
        }
        
    }
        


}