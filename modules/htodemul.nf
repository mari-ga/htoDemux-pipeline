process HTODEMUL{
    publishDir params.outdir, mode:'copy'
    input:
        file preprocess_object
        val quantile_hto
        //For HTODemux
        val kfunc
        val nstarts
        val nsamples
        val nameOutputFileHTO

    output:
        file 'resultHTO.rds'
        file 'resultHTO.csv'

    script:

        """
            Rscript $baseDir/R/HTODemux-args.R --seuratObjectPath $preprocess_object --quantile $quantile_hto --kfunc  $kfunc --nstarts $nstarts --nsamples $nsamples --nameOutputFileHTO $nameOutputFileHTO
        """


}