
process AMPLICOV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/amplicon_coverage_plot:0.3.3--pyhdfd78af_0':
        'biocontainers/amplicon_coverage_plot:0.3.3--pyhdfd78af_0' }"



    input:

    tuple val(meta), path(bam)
    tuple val(meta), val(irma_type)
    path bedpe

    output:

//  output of amplicov is not a sorted bam but just coverage files.

    tuple val(meta), path("*amplicon_coverage.txt"), emit: coverage
    tuple val(meta), path("*amplicon_coverage.html"), emit: coverage_html
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def refname = irma_type

    """
    amplicov --bedpe ${bedpe}  --bam ${bam}  --count_primer -o .  -p ${bam.BaseName}  -r ${refname}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amplicov: \$(amplicov --version |& sed '1!d ; s/amplicov //')
    END_VERSIONS

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${meta}/
    touch ${meta}/stub.amplicon_coverage.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amplicov: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
