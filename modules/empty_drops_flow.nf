params.outdir = 'results'

include { EMPTY_DROPS } from './empty_drops'

workflow EMPTY_DROPS_FLOW{
    take:
       rna_raw 
       niters
       empty
       lower
       testAmbient
       alpha
       ignore
       nameOutputEmpty


    main:
        if(params.empty_drops_mode == "TRUE")
        {
           EMPTY_DROPS(rna_raw,niters,empty,lower,testAmbient,alpha,ignore,nameOutputEmpty)
           
        }

}