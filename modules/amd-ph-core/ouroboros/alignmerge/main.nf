// OUROBOROS_ALIGNMERGE — gather the per-chunk SAMs and produce the winnowed SAM, its score, and the
// refined reference for the next round.
//
// Third stage of ALIGNSPLIT -> ALIGNCHUNK -> ALIGNMERGE. winnow_pileup.py STREAMS the chunk SAMs in
// deterministic order (c0000.sam, c0001.sam, …) in ONE pass: it winnows best-per-read, sums the
// ASM_SCORE (the recursion's convergence objective), and builds the per-position pileup JSON. No giant
// intermediate concatenated SAM is written — that used to be a full write plus re-read of the largest
// file in the pipeline (tens of GB on a deep sample). Each read is wholly within one chunk, so
// per-query records stay contiguous across chunk boundaries, which is the winnow precondition.

process OUROBOROS_ALIGNMERGE {
    tag "${meta.id}_${gene}_r${round}"
    label 'process_medium'

    container 'oamd-bio-ouroboros:1.0.0'

    input:
    tuple val(meta), val(gene), val(round), path(ref), path(chunk_sams)
    val aligner
    val ins_fold_freq
    val del_fold_freq
    val mut_min_strand_frac

    output:
    tuple val(meta), val(gene), val(round), path(ref), path("${gene}.sam"), env('ASM_SCORE')             , emit: sam
    tuple val(meta), val(gene), val(round), path("next-${gene}-r${round}.ref"), env('CONVERGED'), env('INDELS_FOLDED'), emit: refined
    tuple val(meta), val(gene), val(round), path("${gene}.json")                                         , emit: pileup
    path "versions.yml"                                                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sf = (aligner == 'RAMMAP') ? '--score-field 14' : ''
    // Mutation guard: a reference-CHANGING base must clear the variant-calling bar, because mutating the
    // reference is tantamount to calling the variant. Rejects one-strand / thin-data mutations.
    def mut_opt = (mut_min_strand_frac as float) > 0 ? "--mut-min-strand-frac ${mut_min_strand_frac}" : ''
    """
    ASM_SCORE=\$(winnow_pileup.py ${sf} ${ref} ${gene}.sam ${gene}.json --in-sams \$(ls c[0-9]*.sam | sort))

    # Refine folded in: amend the reference from the pileup. The round-unique output name keeps the
    # redirect from clobbering the staged input ref, which IS the prior round's output for round >= 2.
    combine_sam_stats.py ${ref} \\
        --name ${gene} \\
        --insertion-threshold ${ins_fold_freq} \\
        --deletion-threshold ${del_fold_freq} \\
        ${mut_opt} \\
        --indel-flag-file indel.flag \\
        ${gene}.json > next-${gene}-r${round}.ref

    if cmp -s ${ref} next-${gene}-r${round}.ref; then CONVERGED=1; else CONVERGED=0; fi
    INDELS_FOLDED=\$(cat indel.flag)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | grep -oE '[0-9.]+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${gene}.sam ${gene}.json next-${gene}-r${round}.ref
    ASM_SCORE=0
    CONVERGED=1
    INDELS_FOLDED=0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | grep -oE '[0-9.]+' | head -1)
    END_VERSIONS
    """
}
