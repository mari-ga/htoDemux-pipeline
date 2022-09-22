process HTO_DEMUX{
    publishDir params.outdir, mode:'copy'
    label "seurat_process"
    
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

        val quantile_hto
        //For HTODemux
        val kfunc
        val nstarts
        val nsamples
        val seed
        val init
        val nameOutputFileHTO

    output:
        file 'resultHTO.rds'
        file 'resultHTO.csv'

    script:

        """
            Rscript $baseDir/R/HTODemux-complete.R --fileUmi $umi_counts --fileHto $hto_matrix --ndelim $ndelim --selectMethod $selection_method --numberFeatures $number_features --assay $assay --assayName $assayName --margin $margin --normalisationMethod $normalisation_method --quantile $quantile_hto --kfunc  $kfunc --nstarts $nstarts --nsamples $nsamples --seed $seed --init $init --nameOutputFileHTO $nameOutputFileHTO
        """


}