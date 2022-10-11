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
    output:
        file 'emptyDropletsHashed.csv'
        
    script:

        """
          Rscript $baseDir/R/empty_drops.R --fileUmi $hto_raw --niters $niters --empty $empty --lower $lower --testAmbient $testAmbient --alpha $alpha_empty --ignore $ignore --nameOutputEmpty $nameOutputEmpty 
        """

}