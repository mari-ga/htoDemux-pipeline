process EMPTY_REPORT{
    publishDir params.outdir, mode:'copy'
    label "general_results"
    input:
        each classif from columns_classification
        each assig from columns_assignment
        each name from report_names
        val number_reports
      

    output:
        file '*.csv'
        
    script:

        """
          python $baseDir/Python/create_empty_report.py --columns_classification $classif --columns_assignment $assig --number_reports $number_reports --output_report $name
        
        """
}