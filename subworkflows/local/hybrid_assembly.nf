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

    // Seperate reads channel by its type into its channels
    hybrid_dataset.multiMap { it ->
        longreads: [ it[0], it[1] ]
        shortreads: [ it[0], it[2] ]
    }.set { reads }

    // Assemble longreads 
    FLYE (
        reads.longreads
    )
    ch_ver = ch_ver.mix(FLYE.out.versions)
    
    // Branch Assemblied fasta to two channels
    FLYE.out.fasta.multiMap { it ->
        minimap2_ref: it
        asm: it
    }.set { flye_out }

    // Map Assemblied fasta channel with shortreads channel by its meta
    minimap2_ch = reads.shortreads.join(flye_out.minimap2_ref)
    
    // Align shortreads with its corresponding assemblied fasta 
    MINIMAP2_ALIGN (
        minimap2_ch.map{ [ it[0], it[1] ] },
        minimap2_ch.map{ it[2] }
    ) 
    ch_ver = ch_ver.mix(MINIMAP2_ALIGN.out.versions)

    // Sort, index, and convert to BAM
    SAMTOOLS (
       MINIMAP2_ALIGN.out.sam 
    )
    ch_ver = ch_ver.mix(SAMTOOLS.out.versions)

    // Map BAM channel with assemblied fasta channel by its meta
    flye_ch = SAMTOOLS.out.bam.join(flye_out.asm)

    // Correct bases of assemblied fasta using it corresponding BAM files
    PILON (
        flye_ch.map{ [ it[0], it[1], it[2] ] },
        flye_ch.map{ it[3] }
    )
    ch_ver = ch_ver.mix(PILON.out.versions)
    
    fasta = PILON.out.fasta
    
    emit:
    fasta
    versions = ch_ver
}
