// OUROBOROS_INFLATE — expand deflated read patterns back to real reads.
//
// Preprocessing collapses identical reads into patterns plus a deflation index; everything upstream of
// assembly works on patterns. This expands them again so downstream alignment sees TRUE per-read depth.
// Getting this wrong is silent: pattern counts understate depth, and any depth reported from the
// pattern-space side of the pipeline is wrong rather than merely approximate.

process OUROBOROS_INFLATE {
    tag "${meta.id}_${gene}"
    label 'process_low'

    container 'oamd-bio-ouroboros:1.0.0'

    input:
    tuple val(meta), val(gene), path(reads), path(deflation_index)

    output:
    tuple val(meta), val(gene), path("${gene}.fastq"), emit: fastq
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat ${reads} > ${gene}_all.fa
    irma-core xflate --inflate ${deflation_index} ${gene}_all.fa > ${gene}.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma_core: \$(irma-core --version 2>&1 | grep -oE '[0-9][^ ]*' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${gene}.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma_core: \$(irma-core --version 2>&1 | grep -oE '[0-9][^ ]*' | head -1)
    END_VERSIONS
    """
}
