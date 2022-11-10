process HTO_VISUALISATION{

    publishDir path: "$projectDir/graphs", mode:'copy'
    label "seurat_process"
    
    input:
    file result_object
    val assayName
    //Ridge plot params
    val ridgePlot
    val ridgeNCol
    //Scatter features params
    val featureScatter
    val scatterFeat1
    val scatterFeat2
    //Violin plot params
    val vlnplot
    val vlnFeatures
    val vlnLog
    //tSNE
    val tsne
    val tseIdents
    val tsneInvert
    val tsneVerbose
    val tsneApprox
    val tsneDimMax
    val tsePerplexity
    //Heatmap
    val heatmap
    val heatmapNcells

    output:
    //multiples files * with png format
    file '*.png'

    script:
    """
        Rscript $baseDir/R/HTODemux-visualisation.R --pbcmHashtagPath $result_object --assayName $assayName --ridgePlot $ridgePlot --ridgeNCol $ridgeNCol --featureScatter $featureScatter --scatterFeat1 $scatterFeat1 --scatterFeat2 $scatterFeat2 --vlnplot $vlnplot --vlnFeatures $vlnFeatures --vlnLog $vlnLog --tsne $tsne --tseIdents $tseIdents --tsneInvert $tsneInvert --tsneVerbose $tsneVerbose --tsneApprox $tsneApprox --tsneDimMax $tsneDimMax --tsePerplexity $tsePerplexity --heatmap $heatmap --heatmapNcells $heatmapNcells
    """


}