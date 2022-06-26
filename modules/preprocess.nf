

process PREPROCESS{
    publishDir params.outdir, mode:'copy'
    
    input:
        path umi_counts
        path hto_matrix
        val selection_method
        val number_features
        val assay
        val assayName
        val margin
        val normalisation_method
        val nameOutputFile
    output:
        file 'object.rds'

    script:
    """
        Rscript $baseDir/R/pre_processing.R --fileUmi $umi_counts --fileHto $hto_matrix --selectMethod $selection_method --numberFeatures $number_features --assay $assay --assayName $assayName --margin $margin --normalisationMethod $normalisation_method --nameOutputFile $nameOutputFile
    """


}