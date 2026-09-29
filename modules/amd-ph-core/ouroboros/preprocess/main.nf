// OUROBOROS_PREPROCESS — QC, optional adapter trim, and deflation.
//
// Deflation collapses identical reads into patterns plus an index, which is what makes the recruit
// loop affordable: the loop then aligns patterns rather than reads. Everything downstream of gather
// must INFLATE before reporting depth, because a pattern count understates true depth.
//
// No read-level quality trimming: reads are filtered whole, on mean or median quality. Trimming reads
// back would change their alignment footprint and break the parity baseline.

process OUROBOROS_PREPROCESS {
    tag "$meta.id"
    label 'process_medium'

    container 'oamd-bio-python:3.12-orb'

    input:
    // min_len rides in the tuple rather than arriving as a broadcast `val` because it is now
    // derived PER SAMPLE from that sample's own length distribution (see SEQKIT_STATS). A `val`
    // input is one value for the whole run, which is exactly the assumption that put IRMA's
    // 2x150 flu default onto a pre-trimmed library. ph-core forbids custom meta keys, so the
    // per-sample value travels as a tuple element -- the same shape SEQTK_SAMPLE already uses.
    tuple val(meta), path(reads), val(min_len)
    val qual_threshold
    val use_median
    val adapter

    output:
    tuple val(meta), path("*.fa") , emit: reads_fasta
    tuple val(meta), path("*.xfl"), emit: deflation_index
    tuple val(meta), path("QC_log.txt"), emit: qc_log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = reads instanceof List && reads.size() == 2

    def preprocess_opts = "--min-read-quality ${qual_threshold} --min-length ${min_len} --log-file QC_log.txt"
    if (use_median) {
        preprocess_opts += ' --use-median'
    }
    // Adapter trimming is paired-only: the fuzzy match uses the mate to locate the adapter.
    if (paired && adapter) {
        preprocess_opts += " --adapter-trim ${adapter} --a-fuzzy --enforce-clipped-length"
    }

    def read_args = paired ? "${reads[0]} ${reads[1]}" : "${reads}"
    """
    irma-core preprocess \\
        ${preprocess_opts} \\
        ${args} \\
        ${prefix}.xfl \\
        ${read_args} \\
        > ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma_core: \$(irma-core --version | sed 's/irma-core //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa ${prefix}.xfl QC_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma_core: \$(irma-core --version | sed 's/irma-core //')
    END_VERSIONS
    """
}
