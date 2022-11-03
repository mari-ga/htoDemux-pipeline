process HASHED_VISUALISATION{

    publishDir path: "$projectDir/graphs", mode:'copy'
    input:
    file result_object
    val histogram
    val plotLog

    output:
    file '*.png'
    file '*.pdf'

    script:
    """
    Rscript $baseDir/R/hashed-visualisation.R --hashedObjectPath $result_object  --histogram $histogram  --plotLog $plotLog
    """

}