params.outdir = 'results'

include {DEMUXEM_DEMUL} from './demuxem_demul'

workflow DEMUXEM{
    take:
        rna_data
        hto_mat_em
        min_signal
        alpha
        alpha_noise
        tol
        threads
    main:
        if(params.demuxem_mode == "TRUE")
        {
          DEMUXEM_DEMUL(rna_data,hto_mat_em,min_signal,alpha,alpha_noise,tol,threads) 
        }





}