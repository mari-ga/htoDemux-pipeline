process SOLO{
    publishDir params.outdir, mode:'copy'
    label "solo"
    input:
        val rna_file
        val soft
        val max_epochs
        val lr 
        val output_solo 
        
    output:
        path 'solo_prediction.csv', emit: solo_out
        
    script:

        """
          python $baseDir/Python/solo_demul.py --rna_file $rna_data --soft $soft --max_epochs $max_epochs --lr $lr --output $output_sol
        """

}