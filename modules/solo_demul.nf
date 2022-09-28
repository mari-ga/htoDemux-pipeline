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
        if(params.solo_mode == "TRUE")
        {
           SOLO(rna_data,soft,max_epochs,lr,output_solo)
           
        }
    emit:
        SOLO.out
}