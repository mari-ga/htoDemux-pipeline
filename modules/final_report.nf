process FINAL_REPORT{
    publishDir params.outdir, mode:'copy'
    label "general_results"
    input:
        path general_classification
        path solo_prediction
        val output_final
    output:
        file 'final_report.csv'
        
    script:

        """
          python $baseDir/Python/join_reports.py --general_classification $general_classification --solo_results $solo_prediction --output_final $output_final
        
        """
}