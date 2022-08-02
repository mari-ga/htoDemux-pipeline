process DEMUXEM_DEMUL{
    publishDir params.outdir, mode:'copy'
    label "demuxem"
    input:
        path rna_data
        path hto_mat_em
        val threads
        val alpha
        val alpha_noise
        val min_signal
        val tol

    output:
        file 'data.obs'

    script:
    """
        python $baseDir/Python/demuxEM_demul.py --rna_data $rna_data --hto_matrix $hto_mat_em --threads $threads --alpha $alpha --alpha_noise $alpha_noise --min_signal $min_signal --tol $tol
    """

}