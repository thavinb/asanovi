//
// Branch input channel by method of assembly to the designated software.
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

    ch_reads = Channel.empty()
    ch_reads = ch_reads.mix(SPADES.out.scaffolds, FLYE.out.fasta, UNICYCLER.out.scaffolds)

    ch_vers = Channel.empty()
    ch_vers = ch_vers.mix(SPADES.out.versions.first(),FLYE.out.versions.first(),UNICYCLER.out.versions.first())

    emit:
    fasta = ch_reads
    versions = ch_vers
}
