params.outdir = 'results'

include { HASH_SOLO } from './hash_solo'

workflow HASH_SOLO_DEMUL{
    take:
        hto_data 
        priors
        output_file
        output_plot


    main:
        if(params.hash_solo_mode == "TRUE")
        {
           HASH_SOLO(hto_data,priors,output_file,output_plot)
        }
    emit:
        HASH_SOLO.out

}