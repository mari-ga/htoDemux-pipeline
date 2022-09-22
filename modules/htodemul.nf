process HTODEMUL{
    publishDir params.outdir, mode:'copy'
    label "seurat_process"
    input:
        file preprocess_object
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
        file '*.csv'

    script:

        """
            Rscript $baseDir/R/HTODemux-args.R --seuratObjectPath $preprocess_object --quantile $quantile_hto --kfunc  $kfunc --nstarts $nstarts --nsamples $nsamples --seed $seed --init $init --nameOutputFileHTO $nameOutputFileHTO
        """


}