

process PREPROCESS{
    publishDir params.outdir, mode:'copy'
    label "seurat_process"
    input:
        val rdsObject
        val umi_counts
        val hto_matrix
        val ndelim
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
        Rscript $baseDir/R/pre_processing.R  --rdsObject $rdsObject --fileUmi $umi_counts --fileHto $hto_matrix --ndelim $ndelim --selectMethod $selection_method --numberFeatures $number_features --assay $assay --assayName $assayName --margin $margin --normalisationMethod $normalisation_method --nameOutputFile $nameOutputFile
    """


}