process EMPTY_DROPS{
    publishDir params.outdir, mode:'copy'
    label "drops"
    input:
        path hto_raw 
        val niters
        val empty
        val lower
        val testAmbient
        val alpha_empty
        val ignore
        val nameOutputEmpty
        val nameObjectEmpty
    output:
        file 'emptyDropletsHashed.csv'
        path 'emptyDropletsObject.rds', emit: empty_drops_object
        
    script:

        """
          Rscript $baseDir/R/empty_drops.R --fileHto $hto_raw --niters $niters --empty $empty --lower $lower --testAmbient $testAmbient --alpha $alpha_empty --ignore $ignore --nameOutputEmpty $nameOutputEmpty --nameObjectEmpty $nameObjectEmpty
        """

}