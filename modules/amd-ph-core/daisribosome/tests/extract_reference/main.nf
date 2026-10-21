// Copies dais-ribosome's own bundled reference set out of the image
// so tests exercise the wrapper against real biological input
// without committing sequence files.
process EXTRACT_REFERENCE {
    tag "${meta.id}"
    label 'process_single'

    container 'ghcr.io/cdcgov/dais-ribosome:v2.1.0'

    input:
    tuple val(meta), val(annotation_module)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp /app/ribosome_res/${annotation_module}/${annotation_module}-references.fasta ${prefix}.fasta
    """

    stub:
    // Named distinctly from the script block's prefix: nf-test's setup-dependency
    // loader (unlike a normal run) puts both branches in one scope.
    def stubPrefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${stubPrefix}.fasta
    """
}

// Removes the ctype annotation from each FASTA header while leaving sequence
// lines unchanged, allowing tests to exercise --assume-default-ctype.
process EXTRACT_UNANNOTATED_REFERENCE {
    tag "${meta.id}"
    label 'process_single'

    container 'ghcr.io/cdcgov/dais-ribosome:v2.1.0'

    input:
    tuple val(meta), val(annotation_module)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cut -d '|' -f 1 \
        /app/ribosome_res/${annotation_module}/${annotation_module}-references.fasta \
        > ${prefix}.fasta
    """

    stub:
    def stubPrefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${stubPrefix}.fasta
    """
}
