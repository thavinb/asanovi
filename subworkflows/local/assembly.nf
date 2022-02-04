//
// Distribute input by the method to different assembly programs
// Shortread -> SPADES
// Longread -> FLYE
// Hybrid -> UNICYCLER
// 

params.options = [:]

include { SPADES } from '../../modules/nf-core/modules/spades/main.nf' addParams( options: params.options )
include { FLYE } from '../../modules/local/flye' addParams( options: params.options )
include { UNICYCLER } from '../../modules/nf-core/modules/unicycler/main' addParams( options: params.options )

workflow ASSEMBLY {
    take: 
    reads //[meta, [fastq], [fastq_1, fastq_2]]

    main:
    reads.branch { 
        shortread:  it[0].method == "shortread"
        longread:   it[0].method == "longread"
        hybrid:     it[0].method == "hybrid"
    }.set { method }

    SPADES {
        method.shortread
    }
    FLYE {
        method.longread 
    }
    UNICYCLER {
        method.hybrid
    }
    

    emit:
    /* SPADES.out */
    /* FLYE.out */
    /* UNICYCLER.out */
    [SPADES.out, FLYE.out, UNICYCLER.out]
}
