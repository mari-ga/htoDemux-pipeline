process ASSIGNMENT{
    publishDir params.outdir, mode:'copy'
    label "general_results"
    input:
        path hto_assignment
        path multi_assignment
        path hash_drops_assignment
        path hash_solo_assignment
        path demuxem_assignment
        val output_assignment
    output:
        file 'general_assignment.csv'
        
    //
    // 
    script:

        """
          python $baseDir/Python/general_assignment.py --htodemul_assignment $hto_assignment --multiseq_results $multi_assignment  --hashed_drop_results $hash_drops_assignment --hashsolo_results $hash_solo_assignment --demuxem_results $demuxem_assignment --output_assignment $output_assignment
        
        """
}