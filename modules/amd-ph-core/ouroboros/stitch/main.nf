// OUROBOROS_STITCH — place the deep assembled consensus into the full reference coordinate frame.
//
// Iterative assembly produces high-quality calls but TRIMS uncovered columns, so its consensus is
// shorter than the reference and carries no genome coordinates. This projects those calls back onto the
// full panel record, distinguishing three cases that a naive projection would conflate:
//   base   where the assembly aligned one
//   N      where the assembly never spanned, OR an assembly gap whose read depth is below min_depth
//          (a dropout: we do not know what is there)
//   '-'    where an assembly gap IS covered (a real deletion: we know something is missing)
//
// That N-versus-dash distinction is the point. Reuses the reads-on-full-reference BAM and record that
// OUROBOROS_AMPLICONCALL already built rather than remapping.

process OUROBOROS_STITCH {
    tag "${meta.id}_${gene}"
    label 'process_low'

    container 'oamd-bio-python:3.12-orb'

    input:
    tuple val(meta), val(gene), path(amended), path(bam), path(bai), path(ref)
    val min_depth

    output:
    tuple val(meta), val(gene), path("${gene}.stitched.fa"), emit: stitched
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    stitch_consensus.py ${ref} ${amended} \\
        --bam ${bam} \\
        --min-depth ${min_depth} \\
        -N ${gene} \\
        -o ${gene}.stitched.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python3 -c 'import pysam; print(pysam.__version__)' 2>/dev/null || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    touch ${gene}.stitched.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python3 -c 'import pysam; print(pysam.__version__)' 2>/dev/null || echo unknown)
    END_VERSIONS
    """
}
