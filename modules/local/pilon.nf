process PILON {
    tag "$meta.id"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::pilon=1.24" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pilon:1.24--hdfd78af_0':
        'quay.io/biocontainers/pilon:1.24--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem = task.ext.memory ? task.ext.memory.toString().find( /\d+/ ): task.memory.toString().find( /\d+/ )
    """
    export _JAVA_OPTIONS="-Xms${mem}g -Xmx${mem}g"
    pilon \\
        $args \\
        --genome ${reference} \\
        --frags ${bam} \\
        --output ${prefix}.assembly.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pilon: \$(echo \$(pilon --version 2>&1) | grep -e '[0-9]\\.[0-9]*' -o)
    END_VERSIONS
    """
}
