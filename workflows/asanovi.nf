/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowAsanovi.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
/* def modules = params.modules.clone() */

//
// MODULE: Local to the pipeline
//

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' 
include { HYBRID_ASSEMBLY } from '../subworkflows/local/hybrid_assembly.nf' 

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { SPADES as SHORTREAD_ASSEMBLY_SPADES } from '../modules/nf-core/modules/spades/main.nf'
include { FLYE as LONGREAD_ASSEMBLY_FLYE} from '../modules/local/flye' 

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ASANOVI {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // Branching input channel by meta.method 
    //

    INPUT_CHECK.out.reads
    .branch { 
        shortread:  it[0].method == "shortread"
        longread:   it[0].method == "longread"
        hybrid:     it[0].method == "hybrid"
    }.set { types }

    //
    // SHORTREAD ASSEMBLY
    //

    SHORTREAD_ASSEMBLY_SPADES {
        types.shortread
    }    
    ch_versions = ch_versions.mix(SHORTREAD_ASSEMBLY_SPADES.out.versions)

    //
    // LONGREAD ASSEMBLY
    //

    LONGREAD_ASSEMBLY_FLYE {
        types.longread
    }
    ch_versions = ch_versions.mix(LONGREAD_ASSEMBLY_FLYE.out.versions)
    
    //
    // HYBRID ASSEMBLY
    //

    HYBRID_ASSEMBLY (
        types.hybrid
    )
    ch_versions = ch_versions.mix(HYBRID_ASSEMBLY.out.versions)
    
    //
    // MODULE: Pipeline reporting
    //

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')

    )

}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
