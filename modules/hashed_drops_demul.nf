process HASHED_DROPS_DEMUL{
    publishDir params.outdir, mode:'copy'
    
    input:
        
        val hto_matrix
        val nameOutputFileDrops
        val nameOutputFileHashed
        val ambient
        val minProp
        val pseudoCount
        val constAmbient
        val doubletNmads
        val doubletMin
        val confidenMin
        val confidentNmads
        val combinations
    

    output:
        file 'resultHashed_object.rds'
        path 'resultHashed.csv', emit: hashed_drops_results
     

    script:
        """ 
            Rscript $baseDir/R/dropletUtils.R  --fileHto $hto_matrix --nameOutputFileDrops $nameOutputFileDrops --nameOutputFileHashed $nameOutputFileHashed --ambient $ambient --minProp $minProp --pseudoCount $pseudoCount --constAmbient $constAmbient --doubletNmads $doubletNmads --doubletMin $doubletMin --confidenMin $confidenMin --confidentNmads $confidentNmads --combinations $combinations 
        """

}