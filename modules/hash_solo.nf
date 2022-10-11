process HASH_SOLO{
    publishDir params.outdir, mode:'copy'
    label "hash_solo"
    input:
        path hto_data 
        val priors_negative
        val priors_singlet
        val priors_doublet
        val output_file
        val output_plot
    output:
        path 'hash_solo.csv', emit: hash_solo_results
        file 'hash_solo_plot.jpg'

    script:

        """
          python $baseDir/Python/hash_solo.py --hto_data $hto_data --priors_negative $priors_negative --priors_singlet $priors_singlet --priors_doublet $priors_doublet  --output_file $output_file --output_plot $output_plot
        """

}