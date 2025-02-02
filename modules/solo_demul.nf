params.outdir = 'results'

include { SOLO } from './solo'

workflow SOLO_DEMUL{
    take:
        
        rna_data 
        soft
        max_epochs
        lr 
        output_solo 


    main:
           SOLO(rna_data,soft,max_epochs,lr,output_solo)
           solo_ch = SOLO.out.solo_out
        
    emit:
        SOLO_OUT = solo_ch
}