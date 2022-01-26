//
// Check input samplesheet and get read channels
//

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channels(it) }
        .set { reads }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id         = row.sample
    meta.platform   = row.platform

    def array = []
    // TODO Check if file exist for each condition
    if (!file(row.fastq).exists() && !file(row.fastq_1).exists() && !file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }

    if (meta.platform == "shortread") {
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    } else if (meta.platform == "longread") {
        array = [ meta, [ file(row.fastq )] ]
    } else if (meta.platform == "hybrid") {
        array = [ meta, [ file(row.fastq), file(row.fastq_1), file(row.fastq_2) ] ]
    }

    else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return array
}
