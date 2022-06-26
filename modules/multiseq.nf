process MULTISEQ{
    publishDir params.outdir, mode:'copy'
    input:
        file preprocess_object
        val quantile_multi
        //For Multi-seq
        val autoThresh
        val maxiter
        val qrangeFrom
        val qrangeTo
        val qrangeBy
        val verbose
        val nameOutputFileMulti

    output:
        file 'resultMulti.rds'
        file 'resultMulti.csv'

    script:
    
    """
    Rscript $baseDir/R/MULTI-seq.R --seuratObjectPath $preprocess_object --quantile $quantile_multi --autoThresh $autoThresh --maxiter $maxiter --qrangeFrom $qrangeFrom --qrangeTo $qrangeTo --qrangeBy $qrangeBy --verbose $verbose  --nameOutputFileMulti $nameOutputFileMulti
    """


}