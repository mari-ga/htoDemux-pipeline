process EMPTY_DROPS{
    publishDir params.outdir, mode:'copy'
    label "drops"
    input:
        path rna_raw 
        val niters
        val empty
        val lower
        val testAmbient
        val alpha
        val ignore
        val nameOutputEmpty
    output:
        file 'emptyDropletsHashed.csv'
        
    script:

        """
          Rscript $baseDir/R/empty_drops.R --fileUmi $rna_raw --niters $niters --empty $empty --lower $lower --testAmbient $testAmbient --alpha $alpha --ignore $ignore --nameOutputEmpty $nameOutputEmpty 
        """

}