//
// Branch input channel by method of assembly to the designated software.
// 

include { FLYE } from '../../modules/local/flye'
include { MINIMAP2_ALIGN  } from '../../modules/nf-core/modules/minimap2/align/main'
include { SAMTOOLS_SORT as SAMTOOLS } from '../../modules/nf-core/modules/samtools/sort/main'
include { PILON } from '../../modules/local/pilon'

workflow HYBRID_ASSEMBLY{
    take: 
    hybrid_dataset

    main: 
    ch_ver = Channel.empty()

    hybrid_dataset.multiMap { it ->
        longreads: [ it[0], it[1] ]
        shortreads: [ it[0], it[2] ]
    }.set { reads }

    FLYE (
        reads.longreads
    )
    ch_ver = ch_ver.mix(FLYE.out.versions)

    MINIMAP2_ALIGN (
        reads.shortreads,
        FLYE.out.fasta.map{ it[1] }
    ) 
    ch_ver = ch_ver.mix(MINIMAP2_ALIGN.out.versions)

    SAMTOOLS (
       MINIMAP2_ALIGN.out.sam 
    )
    ch_ver = ch_ver.mix(SAMTOOLS.out.versions)

    PILON (
        SAMTOOLS.out.bam,
        FLYE.out.fasta.map{ it[1] }
    )
    ch_ver = ch_ver.mix(PILON.out.versions)
    
    fasta = PILON.out.fasta
    
    emit:
    fasta
    versions = ch_ver
}
