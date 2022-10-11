process REPORT_ASSIGNMENT{
    publishDir params.outdir, mode:'copy'
    label "general_results"
    input:
        val report_name
        
    output:
        file 'general_assignment.csv'
        
    script:

        """
          python $baseDir/Python/create_report.py --report_name $report_name 
        
        """
}