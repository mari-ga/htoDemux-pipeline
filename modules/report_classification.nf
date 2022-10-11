process REPORT_CLASSIFICATION{
    publishDir params.outdir, mode:'copy'
    label "general_results"
    input:
        val report_name
        
    output:
        file 'general_classification.csv'
        
    script:

        """
          python $baseDir/Python/create_report.py --report_name $report_name 
        
        """
}