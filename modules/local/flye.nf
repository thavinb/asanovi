process FLYE {
    tag "$meta.id"
    label 'process_high'
    /* publishDir "${params.outdir}/flye/${meta.id}", */
    /*     mode: params.publish_dir_mode */
    
    conda (params.enable_conda ? "bioconda::flye=2.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flye:2.9--py38h69e0bdc_0':
        'quay.io/biocontainers/flye:2.9--py38h69e0bdc_0' }"

    input:
    tuple val(meta), path(longreads)

    output:
    tuple val(meta), path("*.fasta.gz"), optional:true,  emit: fasta
    tuple val(meta), path("*.gfa.gz")  , optional:true,  emit: gfa
    tuple val(meta), path("*.gv")      , optional:true,  emit: gv
    tuple val(meta), path("*.txt")     , emit: info
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def nano = meta.platform == "nano" ? "--nano-raw $longreads" : "" 
    def pacbio = meta.platform == "pacbio" ? "--pacbio-raw $longreads" : ""
        
    """
    flye \\
        -t $task.cpus \\
        $args \\
        $nano \\
        $pacbio \\
        --out-dir ./

    mv assembly.fasta ${prefix}.assembly.fasta
    gzip -n ${prefix}.assembly.fasta
    mv assembly_graph.gfa ${prefix}.assembly_graph.gfa
    gzip -n ${prefix}.assembly_graph.gfa
    mv assembly_graph.gv ${prefix}.assembly_graph.gv
    mv assembly_info.txt ${prefix}.assembly_info.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$(echo \$(flye --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
}
