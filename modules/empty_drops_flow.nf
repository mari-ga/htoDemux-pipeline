params.outdir = 'results'

include { EMPTY_DROPS } from './empty_drops'

workflow EMPTY_DROPS_FLOW{
    take:
       hto_raw 
       niters
       empty
       lower
       testAmbient
       alpha_empty
       ignore
       nameOutputEmpty
       nameObjectEmpty


    main:
       
           EMPTY_DROPS(hto_raw,niters,empty,lower,testAmbient,alpha_empty,ignore,nameOutputEmpty,nameObjectEmpty)
           
           empty_ch = EMPTY_DROPS.out.empty_drops_object

           
        
    emit:
        empty_out = empty_ch
}