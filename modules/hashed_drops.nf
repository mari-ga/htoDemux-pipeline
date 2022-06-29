params.outdir = 'results'

include { HASHED_DROPS_DEMUL } from './hashed_drops_demul'

workflow HASHED_DROPS{
    take:
        umi_matrix
        hto_matrix
        nameOutputFileDrops
        nameOutputFileHashed
        ambient
        minProp
        pseudoCount
        constAmbient
        doubletNmads
        doubletMin
        confidenMin
        confidentNmads
        combinations



    main:
        if(params.hashedMode == "TRUE")
        {
           HASHED_DROPS_DEMUL(umi_matrix,hto_matrix,nameOutputFileDrops,nameOutputFileHashed,ambient,minProp,pseudoCount,constAmbient,doubletNmads,doubletMin,confidenMin,confidentNmads,combinations) 
        }





}