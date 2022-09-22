params.outdir = 'results'

include { HASHED_DROPS_DEMUL } from './hashed_drops_demul'
include { HASHED_VISUALISATION } from './hashed_drops_visualisation'

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

        histogram
        plotLog



    main:
        if(params.hashedMode == "TRUE")
        {
           HASHED_DROPS_DEMUL(umi_matrix,hto_matrix,nameOutputFileDrops,nameOutputFileHashed,ambient,minProp,pseudoCount,constAmbient,doubletNmads,doubletMin,confidenMin,confidentNmads,combinations)
           if(params.hashedVisualisation == "TRUE")
           {
            HASHED_VISUALISATION(HASHED_DROPS_DEMUL.out[0],histogram,plotLog)
           }

        }





}