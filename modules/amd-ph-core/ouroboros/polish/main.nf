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
//
// Aligner choice (ADR-0015)
// -------------------------
// RAMMAP is the default. A seven-aligner sweep found six of seven produce the same consensus within
// 1-3 bases at EVERY depth: the plurality-refinement loop converges to the same fixpoint regardless of
// how reads were placed, so the loop -- not the aligner -- determines the answer. rammap ties for
// fastest at 10.5x over local_sw with re-indexing included.
//
// LOCALSW is kept for one specific job: the sub-5% minority-variant layer, where the two genuinely
// differ (300 vs 103 calls at full depth on one bin, every difference below 5% frequency). That is a
// use-case switch, not something the pipeline can infer.
//
// bwa, bwa-mem2, minibwa and bowtie2 were evaluated and deliberately NOT adopted -- they match rammap's
// consensus and are slower or equal, and bowtie2 --end-to-end is actively worse (21 wrong bases where
// every other aligner has 0, because it maps 33% of reads where --local maps 96%). The evaluation
// harness at ~/wip/ouroboros/aligner-eval retains them for re-testing on a new assay.

process OUROBOROS_POLISH {
    tag "${meta.id}_${gene}"
    label 'process_high'

    container 'oamd-bio-python:3.12-orb'

    input:
    tuple val(meta), val(gene), path(reads), path(ref), path(deflation)
    val align_opts
    val max_iter

    output:
    tuple val(meta), val(gene), path("polished.ref"), path("polished.sam"), env('CONVERGE_REASON'), emit: aligned
    tuple val(meta), val(gene), env('INDELS_FOLDED')                      , emit: indels_folded
    tuple val(meta), val(gene), path("iter_timings.tsv")                  , emit: timings
    path "versions.yml"                                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def required = ['aligner', 'sw_match', 'sw_mismatch', 'sw_gap_open', 'sw_gap_extend', 'term_ext',
                    'ins_fold_freq', 'del_fold_freq', 'mut_min_strand_frac']
    def missing = required.findAll { k -> !align_opts.containsKey(k) }
    if (missing) {
        error "OUROBOROS_POLISH: align_opts is missing required key(s): ${missing.join(', ')}"
    }

    def o  = align_opts
    // local_sw takes the gap-OPEN penalty exclusive of the first extension; IRMA's convention is
    // inclusive, so the first extension is subtracted here rather than by the caller.
    def go = (o.sw_gap_open as int) - (o.sw_gap_extend as int)
    // 1-based SAM column holding AS:i:, MEASURED per aligner rather than assumed — tag order differs
    // and a wrong column makes win_score() fall back to the CIGAR M-count silently.
    //   bowtie2 (both modes) 12 | minibwa 13 | bwa / bwa-mem2 / rammap 14 | local_sw 12 (default)
    def AS_COL = [ RAMMAP:14 ]          // local_sw's AS is column 12, winnow_pileup's default
    def sf = AS_COL.containsKey(o.aligner) ? "--score-field ${AS_COL[o.aligner]}" : ''
    def known = ['RAMMAP','LOCALSW']
    if (!known.contains(o.aligner)) {
        error "OUROBOROS_POLISH: unknown aligner '${o.aligner}'. Known: ${known.join(', ')}"
    }
    // Both halves of the mutation guard. The depth half was previously never passed, so
    // mut_min_depth stayed 0 and `too_thin` could never fire: a SINGLE read was enough to mutate the
    // reference, and the mutation then propagated into every later iteration.
    def mut_opt = ((o.mut_min_strand_frac as float) > 0 ? "--mut-min-strand-frac ${o.mut_min_strand_frac} " : '') +
                  ((o.get('mut_min_depth', 0) as int) > 0 ? "--mut-min-depth ${o.mut_min_depth} " : '')
    // Optional, defaults off, read with .get() so existing callers are byte-identical. Keeps
    // zero-coverage columns as 'N' instead of dropping them, so the polished reference stays
    // full-length and every downstream coordinate still refers to the same position.
    // Coerce explicitly: a Nextflow CLI param arrives as the STRING "false", which is truthy in
    // Groovy, so a plain truth test silently enables the flag when the caller asked for it off.
    def maskunc = o.get('mask_uncovered', false).toString().toLowerCase() == 'true' ? '--mask-uncovered' : ''
    // Keep the reference base at uncovered columns so the NEXT iteration can still align there.
    // The reference this loop emits is both the alignment target and the reported reference, and
    // 'N' is only right for the second job: once a window is N nothing anchors in it, so it can
    // never recover. Pair with CALL's min_consensus_depth >= 1, which is what then decides the
    // reported consensus -- alone, this flag reports reference bases with no evidence behind them.
    def carryunc = o.get('carry_uncovered', false).toString().toLowerCase() == 'true' ? '--carry-uncovered' : ''
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

    # Seed the outputs from the input reference so they always exist. A degenerate bin — too few reads,
    # or reads too divergent for the aligner to anchor — may never produce a usable iteration, and that
    # is a real result to report, not a reason to fail the run and take the whole cohort with it.
    cp "\$REF" polished.ref
    : > polished.sam
    printf 'iter\\tseconds\\tindex_s\\talign_s\\tscore\\tconverged\\treason\\n' > iter_timings.tsv

    for iter in \$(seq 1 ${max_iter}); do
        t0=\$(date +%s)
        DONE=0
        REASON=not_converged

        # ── ALIGN (threaded over this bin's reads; the SW kernel is itself SIMD) ────────────────────
        # The index-based aligners must REBUILD their index every iteration, because the reference is
        # what the loop is changing. That cost is timed separately so the "loop with reindexing"
        # question has a number rather than an estimate.
        ti=\$(date +%s)
        case "${o.aligner}" in
          LOCALSW)
            local_sw "\$REF" reads.fastq ${o.sw_match} ${o.sw_mismatch} ${go} ${o.sw_gap_extend} \\
                ${task.cpus} ${o.term_ext} > raw_\${iter}.sam ;;
          RAMMAP)
            rammap -a -x sr --secondary no -t ${task.cpus} "\$REF" reads.fastq > raw_\${iter}.sam 2>/dev/null ;;
        esac
        talign=\$(date +%s)

        # ── WINNOW + PILEUP + SCORE (one streaming pass, no concatenated intermediate SAM) ──────────
        SCORE=\$(winnow_pileup.py ${sf} "\$REF" win_\${iter}.sam pileup_\${iter}.json --in-sams raw_\${iter}.sam)
        rm -f raw_\${iter}.sam

        # ── REFINE: amend the reference from this iteration's pileup ────────────────────────────────
        combine_sam_stats.py "\$REF" \\
            --name ${gene} \\
            --insertion-threshold ${o.ins_fold_freq} \\
            --deletion-threshold ${o.del_fold_freq} \\
            ${mut_opt} ${maskunc} ${carryunc} \\
            --indel-flag-file indel.flag \\
            pileup_\${iter}.json > next_\${iter}.ref

        # folded is OR-accumulated: ANY iteration folding an indel marks the bin for the stitch gate.
        [ "\$(cat indel.flag)" != "0" ] && INDELS_FOLDED=1

        # GUARD: if this iteration aligned nothing, the amended reference comes back empty. Feeding that
        # to the next iteration's aligner is what actually fails ("no reference found"), so stop here and
        # keep the last good reference. Seen for real on a 68-read bin where the seeded aligner returned
        # every read unmapped while exhaustive SW aligned them — so this is a live path, not a theoretical one.
        if ! grep -q '[ACGTNacgtn]' next_\${iter}.ref 2>/dev/null; then
            DONE=1
            REASON=no_alignments
            t1=\$(date +%s)
            printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "\$iter" "\$(( t1 - t0 ))" \\
                "\$(( ti - t0 ))" "\$(( talign - ti ))" "\$SCORE" "\$DONE" "\$REASON" >> iter_timings.tsv
            break
        fi

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
        printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "\$iter" "\$(( t1 - t0 ))" \\
            "\$(( ti - t0 ))" "\$(( talign - ti ))" "\$SCORE" "\$DONE" "\$REASON" >> iter_timings.tsv

        PREV_SCORE=\$SCORE
        REF=next_\${iter}.ref
        [ "\$DONE" -eq 1 ] && break
        rm -f win_\$(( iter - 1 )).sam 2>/dev/null || true
    done

    # Terminal state, exported so the caller can route a bin that produced nothing away from variant
    # calling rather than feeding it a zero-byte SAM.
    CONVERGE_REASON=\$REASON

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
    CONVERGE_REASON=reference_stable
    printf 'iter\\tseconds\\tindex_s\\talign_s\\tscore\\tconverged\\treason\\n1\\t0\\t0\\t0\\t0\\t1\\treference_stable\\n' > iter_timings.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma_core: \$(irma-core --version 2>&1 | grep -oE '[0-9][^ ]*' | head -1)
        rammap: \$(rammap --version 2>/dev/null | grep -oE '[0-9][^ ]*' | head -1 || echo unknown)
    END_VERSIONS
    """
}
