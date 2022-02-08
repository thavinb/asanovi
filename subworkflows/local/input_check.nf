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
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channels(it) }
        .set { reads }

    emit:
    // channel: [ val(meta), [ reads ] ]
    reads
    versions = SAMPLESHEET_CHECK.out.versions
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id       = row.sample
    meta.method   = row.method
    meta.platform = row.platform

    def files = [:]
    if (row.fastq) { 
        files.longread = [file(row.fastq)]
    }
    if (row.fastq_1 && row.fastq_2) {
        files.shortread = [file(row.fastq_1), file(row.fastq_2)]
    }

    def array = []

    // Check for file existence

    if (meta.method == "shortread") {
        if (!file(row.fastq_1).exists() && !file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> FastQ file does not exist!\n${row.sample}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    } else if (meta.method == "longread") {
        if (!file(row.fastq).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> FastQ file does not exist!\n${row.sample}"
        }
        array = [ meta, [ file(row.fastq) ] ]
    } else if (meta.method == "hybrid") {
        if (!file(row.fastq).exists() && !file(row.fastq_1).exists() && !file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> FastQ file does not exist!\n${row.sample}"
        }
        array = [ meta, files.longread, files.shortread ]
    }

    else {
        exit 1, "ERROR: Please check input samplesheet -> Platform muse be specified!\n${row.sample}" 
    }
    return array
}
