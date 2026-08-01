// OUROBOROS_POLISH — the whole iterative polish loop for ONE bin, in ONE task.
//
// Each iteration aligns the bin's reads to the current reference, winnows to best-per-read, scores the
// alignment, and amends the reference from the resulting pileup. The loop stops when the reference
// reaches a byte fixpoint, when the alignment score plateaus, or at max_iter. The HIGHEST-scoring
// iteration's (reference, SAM) pair is what variant calling consumes — picking the best rather than
// the last handles both a plateau and an outright score decrease without threading a "previous best"
// through any state.
//
// Why the loop is in here
// -----------------------
// The same reason as OUROBOROS_RECRUIT: no in-workflow iteration survives more than one item in
// flight. `until` closes the whole channel, `filter` deadlocks, and `.recurse()` is value-channel only
// — which is why the previous design bundled every bin of a sample into ONE recursion state and then
// needed a skip-converged workaround to stop converged bins re-aligning against an identical reference
// every round.
//
// Looping here makes each BIN an ordinary Nextflow task. Bins run concurrently and converge
// independently, with no round barrier: a shallow contaminant stops at its own convergence instead of
// waiting on the deepest bin every iteration.
//
// What this gives up, and why that is the right trade
// ---------------------------------------------------
// Alignment now uses the aligner's own threads inside the task rather than fanning chunks out across
// Nextflow tasks. Measured, that fan-out fired on 1 of 33 cohort samples and was worth 3.1% of cohort
// wall, while the sequential rounds it replaces are ~36%. Both aligners thread internally (local_sw
// rayon, rammap -t) and their kernels are SIMD, so the parallelism moves rather than disappears —
// from chunks-across-tasks to bins-across-tasks plus threads within.

process OUROBOROS_POLISH {
    tag "${meta.id}_${gene}"
    label 'process_high'

    container 'oamd-bio-ouroboros:1.0.0'

    input:
    tuple val(meta), val(gene), path(reads), path(ref), path(deflation)
    val align_opts
    val max_iter

    output:
    tuple val(meta), val(gene), path("polished.ref"), path("polished.sam"), emit: aligned
    tuple val(meta), val(gene), env('INDELS_FOLDED')                      , emit: indels_folded
    tuple val(meta), val(gene), path("iter_timings.tsv")                  , emit: timings
    path "versions.yml"                                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def required = ['aligner', 'sw_match', 'sw_mismatch', 'sw_gap_open', 'sw_gap_extend', 'term_ext',
                    'ins_fold_freq', 'del_fold_freq', 'mut_min_strand_frac']
    def missing = required.findAll { !align_opts.containsKey(it) }
    if (missing) {
        error "OUROBOROS_POLISH: align_opts is missing required key(s): ${missing.join(', ')}"
    }

    def o  = align_opts
    // local_sw takes the gap-OPEN penalty exclusive of the first extension; IRMA's convention is
    // inclusive, so the first extension is subtracted here rather than by the caller.
    def go = (o.sw_gap_open as int) - (o.sw_gap_extend as int)
    def sf = (o.aligner == 'RAMMAP') ? '--score-field 14' : ''
    def mut_opt = (o.mut_min_strand_frac as float) > 0 ? "--mut-min-strand-frac ${o.mut_min_strand_frac}" : ''
    """
    set -e

    # Inflate ONCE. The reads are iteration-invariant — only the reference moves — so inflating inside
    # the loop would repeat the same expansion every round. Depth from here on is TRUE per-read depth.
    cat ${reads} > all_patterns.fa
    irma-core xflate --inflate ${deflation} all_patterns.fa > reads.fastq

    REF=${ref}
    BEST_SCORE=""
    INDELS_FOLDED=0
    PREV_SCORE=""
    printf 'iter\\tseconds\\tscore\\tconverged\\treason\\n' > iter_timings.tsv

    for iter in \$(seq 1 ${max_iter}); do
        t0=\$(date +%s)
        DONE=0
        REASON=not_converged

        # ── ALIGN (threaded over this bin's reads; the SW kernel is itself SIMD) ────────────────────
        if [ "${o.aligner}" = "RAMMAP" ]; then
            rammap -a -x sr --secondary no -t ${task.cpus} "\$REF" reads.fastq > raw_\${iter}.sam 2>/dev/null
        else
            local_sw "\$REF" reads.fastq ${o.sw_match} ${o.sw_mismatch} ${go} ${o.sw_gap_extend} \\
                ${task.cpus} ${o.term_ext} > raw_\${iter}.sam
        fi

        # ── WINNOW + PILEUP + SCORE (one streaming pass, no concatenated intermediate SAM) ──────────
        SCORE=\$(winnow_pileup.py ${sf} "\$REF" win_\${iter}.sam pileup_\${iter}.json --in-sams raw_\${iter}.sam)
        rm -f raw_\${iter}.sam

        # ── REFINE: amend the reference from this iteration's pileup ────────────────────────────────
        combine_sam_stats.py "\$REF" \\
            --name ${gene} \\
            --insertion-threshold ${o.ins_fold_freq} \\
            --deletion-threshold ${o.del_fold_freq} \\
            ${mut_opt} \\
            --indel-flag-file indel.flag \\
            pileup_\${iter}.json > next_\${iter}.ref

        # folded is OR-accumulated: ANY iteration folding an indel marks the bin for the stitch gate.
        [ "\$(cat indel.flag)" != "0" ] && INDELS_FOLDED=1

        # ── track the BEST-scoring iteration: its (ref, sam) is what CALL consumes ──────────────────
        if [ -z "\$BEST_SCORE" ] || [ "\$SCORE" -gt "\$BEST_SCORE" ]; then
            BEST_SCORE=\$SCORE
            cp "\$REF" polished.ref
            cp win_\${iter}.sam polished.sam
        fi

        # ── CONVERGENCE: byte fixpoint, or a score plateau ──────────────────────────────────────────
        # The plateau test matters for samples whose reference jitters between two states without ever
        # byte-fixpointing; without it those bins would run to max_iter every time.
        if cmp -s "\$REF" next_\${iter}.ref; then
            DONE=1; REASON=reference_stable
        elif [ -n "\$PREV_SCORE" ] && [ "\$SCORE" -le "\$PREV_SCORE" ]; then
            DONE=1; REASON=score_plateau
        fi

        t1=\$(date +%s)
        printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "\$iter" "\$(( t1 - t0 ))" "\$SCORE" "\$DONE" "\$REASON" >> iter_timings.tsv

        PREV_SCORE=\$SCORE
        REF=next_\${iter}.ref
        [ "\$DONE" -eq 1 ] && break
        rm -f win_\$(( iter - 1 )).sam 2>/dev/null || true
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma_core: \$(irma-core --version 2>&1 | grep -oE '[0-9][^ ]*' | head -1)
        rammap: \$(rammap --version 2>/dev/null | grep -oE '[0-9][^ ]*' | head -1 || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    touch polished.ref polished.sam
    INDELS_FOLDED=0
    printf 'iter\\tseconds\\tscore\\tconverged\\treason\\n1\\t0\\t0\\t1\\treference_stable\\n' > iter_timings.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma_core: \$(irma-core --version 2>&1 | grep -oE '[0-9][^ ]*' | head -1)
        rammap: \$(rammap --version 2>/dev/null | grep -oE '[0-9][^ ]*' | head -1 || echo unknown)
    END_VERSIONS
    """
}
