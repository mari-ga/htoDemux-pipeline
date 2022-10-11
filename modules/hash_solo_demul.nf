params.outdir = 'results'

include { HASH_SOLO } from './hash_solo'

workflow HASH_SOLO_DEMUL{
    take:
        hto_data 
        priors_negative
        priors_singlet
        priors_doublet
        output_file
        output_plot


    main:
        if(params.hash_solo_mode == "TRUE")
        {
           HASH_SOLO(hto_data,priors_negative,priors_singlet,priors_doublet,output_file,output_plot)
           hash_solo_ch = HASH_SOLO.out.hash_solo_results
        }
    emit:
        HASH_SOLO_OUT = hash_solo_ch

}