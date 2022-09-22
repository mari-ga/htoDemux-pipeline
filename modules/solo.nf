process SOLO{
    publishDir params.outdir, mode:'copy'
    label "solo"
    input:
        val rna_data 
        val soft
        val max_epochs
        val lr 
        val output_solo 
        
    output:
        file 'solo_prediction.csv'
        
    script:

        """
          python $baseDir/Python/solo.py --rna_data $rna_data --soft $soft --max_epochs $max_epochs --lr $lr --output $output_demux
        """

}