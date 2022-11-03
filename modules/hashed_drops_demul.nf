process HASHED_DROPS_DEMUL{
    publishDir params.outdir, mode:'copy'
    
    input:
        
        path empty_drops_result
        val rawData
        val hashtag_data
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
    

    output:
        file 'resultHashed_object.rds'
        path 'resultHashed.csv', emit: hashed_drops_results
     

    script:
        """ 
            Rscript $baseDir/R/dropletUtils.R  --rdsObject $empty_drops_result --rawData $rawData --fileHto $hashtag_data --nameOutputFileDrops $nameOutputFileDrops --nameOutputFileHashed $nameOutputFileHashed --ambient $ambient --minProp $minProp --pseudoCount $pseudoCount --constAmbient $constAmbient --doubletNmads $doubletNmads --doubletMin $doubletMin --confidenMin $confidenMin --confidentNmads $confidentNmads 
        """

}