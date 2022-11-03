params.outdir = 'results'

include { HASHED_DROPS_DEMUL } from './hashed_drops_demul'
include { HASHED_VISUALISATION } from './hashed_drops_visualisation'
include { EMPTY_DROPS } from './empty_drops'

workflow HASHED_DROPS{
    take:
        //For hashed Drops
        empty_drops_result
        rawData
        hashtag_data
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

        //For hashed Drops Visualisation
        histogram
        plotLog

        //For empty drops - one last time
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
        if(params.rawData == "FALSE"){
            HASHED_DROPS_DEMUL(empty_drops_result,rawData,hashtag_data,nameOutputFileDrops,nameOutputFileHashed,ambient,minProp,pseudoCount,constAmbient,doubletNmads,doubletMin,confidenMin,confidentNmads)
            hashed_drops_ch = HASHED_DROPS_DEMUL.out.hashed_drops_results
            if(params.hashedVisualisation == "TRUE")
            {
            HASHED_VISUALISATION(HASHED_DROPS_DEMUL.out[0],histogram,plotLog)
            }
        }else{
            print("Hashed Drops combined with Empty drops")
            EMPTY_DROPS(hto_raw,niters,empty,lower,testAmbient,alpha_empty,ignore,nameOutputEmpty,nameObjectEmpty)
            HASHED_DROPS_DEMUL(EMPTY_DROPS.out[1],rawData,hashtag_data,nameOutputFileDrops,nameOutputFileHashed,ambient,minProp,pseudoCount,constAmbient,doubletNmads,doubletMin,confidenMin,confidentNmads)
            hashed_drops_ch = HASHED_DROPS_DEMUL.out.hashed_drops_results
            if(params.hashedVisualisation == "TRUE")
            {
            HASHED_VISUALISATION(HASHED_DROPS_DEMUL.out[0],histogram,plotLog)
            }
        }
        
    emit:
        HASHED_DROPS_OUT = hashed_drops_ch





}