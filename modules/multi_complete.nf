process MULTI_SEQ{
    publishDir params.outdir, mode:'copy'
    
    input:
    val umi_counts
    val hto_matrix
    val ndelim
    val selection_method
    val number_features
    val assay
    val assayName
    val margin
    val normalisation_method

    val quantile_multi
    val autoThresh
    val maxiter
    val qrangeFrom
    val qrangeTo
    val qrangeBy
    val verbose
    val nameOutputFileMulti
    val nameClassificationFileMulti

    output:
        file 'resultMulti_object.rds'
        path 'resultMulti.csv', emit: classification_multi

    script:
    """
    
    Rscript $baseDir/R/multi-complete.R --fileUmi $umi_counts --fileHto $hto_matrix --ndelim $ndelim --selectMethod $selection_method --numberFeatures $number_features --assay $assay --assayName $assayName --margin $margin --normalisationMethod $normalisation_method --quantile_multi $quantile_multi --autoThresh $autoThresh --maxiter $maxiter --qrangeFrom $qrangeFrom --qrangeTo $qrangeTo --qrangeBy $qrangeBy --verbose $verbose  --nameOutputFileMulti $nameOutputFileMulti --nameClassificationFileMulti $nameClassificationFileMulti
    """
}