process DAISRIBOSOME {
    tag "${meta.id}"
    label 'process_medium'

    // Multi-platform OCI index with native linux/amd64 and linux/arm64 images.
    container 'ghcr.io/cdcgov/dais-ribosome:v2.1.0'

    input:
    tuple val(meta), path(sequences)
    val annotation_module

    output:
    tuple val(meta), path("*.seq.txt")    , emit: product_seq
    tuple val(meta), path("*.ins.txt")    , emit: product_ins
    tuple val(meta), path("*.del.txt")    , emit: product_del
    tuple val(meta), path("*.gen_seq.txt"), emit: genome_seq
    tuple val(meta), path("*.gen_ins.txt"), emit: genome_ins
    tuple val(meta), path("*.gen_del.txt"), emit: genome_del
    tuple val("${task.process}"), val('dais-ribosome'), eval("ribosome --version | cut -d ' ' -f 2"), topic: versions, emit: versions_daisribosome

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ribosome \\
        ${args} \\
        --module ${annotation_module} \\
        --verbose \\
        --threads ${task.cpus} \\
        --output-prefix ${prefix} \\
        ${sequences}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch \\
        ${prefix}.seq.txt \\
        ${prefix}.ins.txt \\
        ${prefix}.del.txt \\
        ${prefix}.gen_seq.txt \\
        ${prefix}.gen_ins.txt \\
        ${prefix}.gen_del.txt
    """
}
