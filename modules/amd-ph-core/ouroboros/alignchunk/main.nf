// OUROBOROS_ALIGNCHUNK — align ONE read chunk to the reference. The scatter unit.
//
// Second stage of ALIGNSPLIT -> ALIGNCHUNK -> ALIGNMERGE. One task per chunk; both aligners thread
// internally over the chunk's reads (local_sw rayon -t, rammap -t), and each read's Smith-Waterman
// kernel is itself SIMD. Memory stays flat (local_sw batches, rammap streams), so a chunk is sized to
// soak up most of a sample as ONE multicore task — ALIGNSPLIT only fans out for super-deep samples.
// On a local executor a few fat tasks fill the cores; on AWS Batch they spread across instances.
// Emits the chunk's SAM (with its own header) for ALIGNMERGE to stream and winnow.

process OUROBOROS_ALIGNCHUNK {
    tag "${meta.id}_${gene}_r${round}_c${chunk_id}"
    label 'process_high'

    container 'oamd-bio-ouroboros:1.0.0'

    input:
    tuple val(meta), val(gene), val(round), val(chunk_id), path(ref), path(chunk_fastq)
    val aligner
    val sw_match
    val sw_mismatch
    val sw_gap_open
    val sw_gap_extend
    val term_ext

    output:
    tuple val(meta), val(gene), val(round), path("c${chunk_id}.sam"), emit: sam
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // local_sw takes the gap-OPEN penalty exclusive of the first extension; IRMA's sw_gap_open is
    // inclusive, so the first extension is subtracted here rather than in the caller.
    def go   = (sw_gap_open as int) - (sw_gap_extend as int)
    def args = task.ext.args ?: ''
    if (aligner == 'RAMMAP')
        """
        rammap -a -x sr --secondary no -t ${task.cpus} ${args} ${ref} ${chunk_fastq} > c${chunk_id}.sam 2>/dev/null

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rammap: \$(rammap --version 2>/dev/null | grep -oE '[0-9][^ ]*' | head -1)
        END_VERSIONS
        """
    else
        """
        local_sw ${ref} ${chunk_fastq} ${sw_match} ${sw_mismatch} ${go} ${sw_gap_extend} ${task.cpus} ${term_ext} > c${chunk_id}.sam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            local_sw: \$(rammap --version 2>/dev/null | grep -oE '[0-9][^ ]*' | head -1)
        END_VERSIONS
        """

    stub:
    """
    touch c${chunk_id}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rammap: \$(rammap --version 2>/dev/null | grep -oE '[0-9][^ ]*' | head -1)
    END_VERSIONS
    """
}
