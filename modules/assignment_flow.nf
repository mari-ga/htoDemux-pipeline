params.outdir = 'results'

include { ASSIGNMENT } from './assignment'

workflow ASSIGNMENT_WORKFLOW{
    take:
        hto_assignment
        multi_assignment
        hash_drops_assignment
        hash_solo_assignment
        demuxem_assignment
        output_assignment

    //
    main:
      
           ASSIGNMENT(hto_assignment,multi_assignment,hash_drops_assignment,hash_solo_assignment,demuxem_assignment,output_assignment)
           
        

}