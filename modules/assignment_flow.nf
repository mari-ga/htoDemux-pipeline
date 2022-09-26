params.outdir = 'results'

include { ASSIGNMENT } from './assignment'

workflow ASSIGNMENT_WORKFLOW{
    take:
        hto_assignment
        multi_assignment
        hash_drops_assignment
        demuxem_assignment
        hash_solo_assignment
        output_assignment


    main:
        if(params.general_mode == "TRUE")
        {
           ASSIGNMENT(hto_assignment,multi_assignment,hash_drops_assignment,demuxem_assignment,hash_solo_assignment,output_assignment)
           
        }

}