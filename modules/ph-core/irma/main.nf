process IRMA {
    tag "$meta.id"
    label 'process_high'

    container 'cdcgov/irma:v1.1.5'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}/*")    , emit: results
    path("${meta.id}/amended_consensus/*.fa"), optional: true, emit: amended_consensus
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def module = "${params.irma_module}"

    """
    IRMA \\
        $module \\
        $reads \\
        ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma: \$(echo \$(IRMA | grep -o "v[0-9][^ ]*" | cut -c 2-))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def module = "${params.irma_module}"

    """
    touch ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma: \$(echo \$(IRMA | grep -o "v[0-9][^ ]*" | cut -c 2-))
    END_VERSIONS
    """
}
