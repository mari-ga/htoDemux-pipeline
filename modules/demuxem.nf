process DEMUXEM{
    publishDir params.outdir, mode:'copy'
    label "demuxem"
    input:
        path rna_data 
        path hto_data
        val alpha
        val alpha_noise
        val tol
        val n_threads
        val min_signal
        val output_demux
    output:
        file 'demuxEm.csv'
        
    script:

        """
          python $baseDir/Python/demuxEM_demul.py --rna_data $rna_data --hto_matrix $hto_data --output $output_demux
        """

}