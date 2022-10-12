process CLASSIFICATION{
    publishDir params.outdir, mode:'copy'
    label "general_results"
    input:
        path hto_classification
        path multi_classification
        path hash_drops_classification
        //path demuxem_classification
        path hash_solo_classification
        val output_classification
    output:
        file 'general_classification.csv'
        
    script:

        """
          python $baseDir/Python/general_classification.py --htodemul_classification $hto_classification --multiseq_results $multi_classification --hashed_drop_results $hash_drops_classification --hashsolo_results $hash_solo_classification --output_assignment $output_classification
        
        """
}