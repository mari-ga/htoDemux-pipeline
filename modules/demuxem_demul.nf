params.outdir = 'results'

include { DEMUXEM } from './demuxem'

workflow DEMUXEM_DEMUL{
    take:
        rna_data
        hto_data
        alpha
        alpha_noise
        tol
        n_threads
        min_signal
        output_demux


    main:
        if(params.demuxem_mode == "TRUE")
        {
           DEMUXEM(rna_data,hto_data,alpha,alpha_noise,tol,n_threads, min_signal,output_demux)
           
        }

}