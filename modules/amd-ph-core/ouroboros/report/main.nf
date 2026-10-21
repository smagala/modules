// OUROBOROS_REPORT — per-gene diagnostic PDFs (coverage and heuristics).
//
// Human-facing artifacts only; nothing downstream consumes them. Kept as its own process so the
// plotting stack stays off the compute-heavy calling path. Plot failures are swallowed deliberately:
// a diagram is never a reason to fail a gene, so both outputs are optional.

process OUROBOROS_REPORT {
    tag "${meta.id}_${gene}"
    label 'process_single'

    container 'oamd-bio-python:3.12-orb'

    input:
    tuple val(meta), val(gene), path(coverage), path(variants), path(alleles), path(pairing)
    val min_variant_qual
    val min_variant_freq
    val min_column_coverage
    val min_confidence

    output:
    tuple val(meta), val(gene), path("${gene}-coverageDiagram.pdf"), emit: coverage_diagram , optional: true
    tuple val(meta), val(gene), path("${gene}-heuristics.pdf")     , emit: heuristic_diagram, optional: true
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export MPLCONFIGDIR="\$PWD"
    n=\$(wc -l < ${variants} 2>/dev/null || echo 0)
    # The full coverage diagram overlays minority variants, which needs both the variant table and the
    # paired stats; fall back to the depth-only plot for single-end or no-variant genes.
    if [ "\$n" -gt 1 ] && [ -s "${pairing}" ]; then
        coverage_diagram.py "${meta.id}" ${gene} ${coverage} ${variants} ${pairing} ${gene}-coverageDiagram.pdf 2>/dev/null || true
    else
        coverage_diagram.py "${meta.id}" ${gene} ${coverage} ${gene}-coverageDiagram.pdf 2>/dev/null || true
    fi
    heuristic_diagram.py ${min_variant_qual} ${min_variant_freq} ${min_column_coverage} ${min_confidence} ${alleles} ${gene}-heuristics.pdf 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)' 2>/dev/null || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)' 2>/dev/null || echo unknown)
    END_VERSIONS
    """
}
