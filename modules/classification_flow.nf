params.outdir = 'results'

include { CLASSIFICATION } from './classification'

workflow CLASSIFICATION_WORKFLOW{
    
    take:
        hto_classification
        multi_classification
        hash_drops_classification
        hash_solo_classification
        demuxem_classification
        output_classification


    main:
      
           CLASSIFICATION(hto_classification,multi_classification,hash_drops_classification,hash_solo_classification,demuxem_classification,output_classification)
            clasification_ch = CLASSIFICATION.out.classification_results
    emit:
        classification_out = clasification_ch

        

}