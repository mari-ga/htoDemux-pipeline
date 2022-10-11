process DEMUXEM{
    publishDir params.outdir, mode:'copy'
    label "demuxem"
    input:
        path rna_raw 
        path hto_matrix
        val alpha
        val alpha_noise
        val tol
        val n_threads
        val min_signal
        val output_demux
    output:
        path 'demuxEm.csv', emit: demuxem_out
        
    script:

        """
          python $baseDir/Python/demuxEM_demul.py --rna_data $rna_raw --hto_matrix $hto_matrix --output $output_demux
        """

}