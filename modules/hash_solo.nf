process HASH_SOLO{
    publishDir params.outdir, mode:'copy'
    label "hash_solo"
    input:
        path hto_data 
        val priors
        val output_file
        val output_plot
    output:
        file 'hash_solo.csv'
        file 'hash_solo.jpg'

    script:

        """
          python $baseDir/Python/hash_solo.py --hto_data $hto_data --output_file $output_file --output_plot $output_plot
        """

}