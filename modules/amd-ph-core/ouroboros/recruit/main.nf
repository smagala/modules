// OUROBOROS_RECRUIT — the whole recruit/refine loop for ONE sample, in ONE task.
//
// Each round is MATCH -> SORT -> per-gene ALIGN+REFINE -> convergence close, and the loop turns inside
// the task until the reference stops moving, nothing is left to recruit, or max_rounds is reached.
//
// Why the loop is in here rather than in the workflow
// ---------------------------------------------------
// Nextflow's in-workflow iteration options all fail once more than one item is in flight:
//
//   * topic feedback + `until` — `until` CLOSES the channel rather than dropping the matching item, so
//     the first sample to converge truncates the loop of every other sample still running. Silently:
//     truncated recruitment yields a short reference, then a bad assembly, with no error anywhere.
//   * topic feedback + `filter` — drops items individually, but the topic never closes, so it deadlocks.
//   * `.recurse()` — value-channel only, which forces every bin of a sample into one bundled state.
//
// The remaining alternative was to serialise each bin's boundary to a directory and dispatch a separate
// `nextflow run` per bundle. That works, but it makes the filesystem the channel and writes the read
// data twice.
//
// Looping in-task avoids all of it: per-sample independence (and, downstream, per-bin) becomes an
// ordinary Nextflow scatter of one task each, and no state touches disk between rounds.
//
// It also lands ADR-0013 from the other direction: across 34 traced samples R1 was 97.3% of gather
// compute, while all 43 tail tasks combined did 25.5 s of work and waited 8-34x longer than they
// worked. Those task boundaries bought nothing.
//
// Observability: per-round recruitment counts and timings are emitted as DATA (sorted_read_stats.txt,
// round_timings.tsv) rather than being recovered from execution_trace rows. That survives any task
// granularity and is analysable across a cohort, which trace rows are not.

process OUROBOROS_RECRUIT {
    tag "${meta.id}"
    label 'process_high'

    container 'oamd-bio-python:3.12-orb'

    input:
    tuple val(meta), path(reads), path(refs)
    val recruit_opts
    val max_rounds

    output:
    tuple val(meta), path("R*-*.ref")              , emit: gene_refs , optional: true
    tuple val(meta), path("R*-*.fa")               , emit: sorted    , optional: true
    tuple val(meta), path("nomatch.fa")            , emit: nomatch
    tuple val(meta), path("chimeric.fa")           , emit: chimeric
    tuple val(meta), path("sorted_read_stats.txt") , emit: sort_stats
    tuple val(meta), path("round_timings.tsv")     , emit: timings
    tuple val(meta), path("R*.defer.tsv")          , emit: deferred  , optional: true
    tuple val(meta), path("R*.collapse.tsv")       , emit: collapse  , optional: true
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Settings arrive as one Map rather than twenty-one positional val inputs, which would mis-bind
    // silently on any reorder. Validate here so a missing key fails with a name instead of
    // interpolating "null" into a shell command.
    def required = ['match_aligner', 'gather_aligner', 'blast_task', 'blast_dust', 'blast_evalue',
                    'blast_min_hsp_len', 'blast_min_hsp_pid', 'blast_incl_chim', 'min_recruit_count',
                    'min_recruit_patterns', 'sort_groups', 'ban_groups', 'secondary_min_depth',
                    'secondary_min_frac', 'min_ca', 'min_fa', 'skip_elongation', 'gather_position_bias',
                    'gather_del_type', 'min_column_coverage', 'gather_ref_floor']
    def missing = required.findAll { k -> !recruit_opts.containsKey(k) }
    if (missing) {
        error "OUROBOROS_RECRUIT: recruit_opts is missing required key(s): ${missing.join(', ')}"
    }

    def o         = recruit_opts
    def incl      = o.blast_incl_chim ? '--incl-chim' : ''
    def sg        = o.sort_groups ? "-P '${o.sort_groups}'" : ''
    def bg        = o.ban_groups  ? "-B '${o.ban_groups}'"  : ''
    // Secondary-target keep (co-infection / contamination): prefer a COVERAGE floor on the secondary's
    // own reference, which is invariant to how deep the primary is, over the fraction-of-primary gate,
    // which moves with the primary's depth. \$REFS and \$RL resolve in the script body.
    def sec       = (o.secondary_min_depth && (o.secondary_min_depth as float) > 0)
                        ? "--ref-fasta \$REFS --read-len \$RL --secondary-min-depth ${o.secondary_min_depth}"
                        : ((o.secondary_min_frac && (o.secondary_min_frac as float) > 0)
                            ? "--secondary-frac ${o.secondary_min_frac} --secondary-floor ${o.min_ca}"
                            : '')
    // Ambiguity deferral is OPTIONAL and defaults off, so it is read with .get() rather than added to
    // the required-key list — existing callers keep byte-identical behaviour without being edited.
    // Only the rammap MATCH path has it; the BLAST path is untouched.
    def defer     = o.get('defer_ambiguous', false)
                        ? "--defer-ambiguous --margin-matches ${o.get('margin_matches', 2)} " +
                          "--margin-frac ${o.get('margin_frac', 0.02)} --defer-report R\${round}.defer.tsv"
                        : ''
    // Optional, defaults off, read with .get() like defer_ambiguous so existing callers are
    // unchanged. \$REFS and the round's hit file resolve in the script body.
    def cgroups   = o.get('collapse_groups', null)
    def collapse  = ''
    if (cgroups) {
        def gl = (cgroups instanceof CharSequence) ? [cgroups as String] : (cgroups as List)
        def hitsFile = o.match_aligner == 'BLAST' ? 'hits.tsv' : 'hits.paf'
        def hitsFmt  = o.match_aligner == 'BLAST' ? 'blast6'    : 'paf'
        collapse = gl.collect { "--collapse-group '${it}'" }.join(' ') +
                   " --hits ${hitsFile} --hits-format ${hitsFmt}" +
                   " --hits-min-len ${o.blast_min_hsp_len} --hits-min-pid ${o.blast_min_hsp_pid}" +
                   " --breadth-min-depth ${o.get('collapse_breadth_min_depth', 1)}" +
                   " --breadth-tol ${o.get('collapse_breadth_tol', 0.02)}" +
                   " --collapse-min-breadth-frac ${o.get('collapse_min_breadth_frac', 0.5)}" +
                   " --collapse-warn-breadth ${o.get('collapse_warn_breadth', 0.50)}" +
                   " --collapse-on-suspect ${o.get('collapse_on_suspect', 'skip')}" +
                   " --collapse-report R\${round}.collapse.tsv"
        // Breadth needs per-target lengths; the secondary gate may already have passed --ref-fasta.
        if (!sec.contains('--ref-fasta')) {
            collapse = "--ref-fasta \$REFS " + collapse
        }
    }
    def elo       = o.skip_elongation ? '--skip-elongation' : ''
    def tw        = o.gather_position_bias ? 3 : 0
    def del_ambig = o.gather_del_type == 'NNN' ? '--delete-by-ambiguity' : ''
    def pb        = o.gather_position_bias ? '--position-bias-frac 0.9 --position-bias-min-span 5' : ''
    def mrd       = (o.gather_ref_floor && o.gather_del_type == 'REF' && o.min_column_coverage > 0)
                        ? "--min-ref-depth ${o.min_column_coverage}" : ''
    """
    set -e

    # Round state lives in this task's own working directory. There is no instate/outstate handoff,
    # because there is no task boundary between rounds: REFS and READS advance each iteration, and the
    # cumulative per-gene JSONs accumulate as R*-<gene>.json without being copied anywhere.
    REFS=${refs}
    READS=${reads}
    LAST_ROUND=1
    printf 'round\\tseconds\\tconverged\\treason\\n' > round_timings.tsv

    for round in \$(seq 1 ${max_rounds}); do
        t0=\$(date +%s)
        LAST_ROUND=\$round
        DONE=0
        REASON=not_converged

        # Round 1 applies the real recruitment floors; later rounds accept any read that matches,
        # because by then the reference has moved and a single read is meaningful evidence.
        if [ "\$round" -eq 1 ]; then
            MRC=${o.min_recruit_count}; MRP=${o.min_recruit_patterns}
        else
            MRC=1; MRP=1
        fi

        # ── MATCH ───────────────────────────────────────────────────────────────────────────────────
        # A round can inherit no carry reads (a shallow run where round 1 recruited everything). rammap
        # spins indefinitely on a zero-byte query rather than returning empty, so the aligner must never
        # see one: skip MATCH and let the empty-class branch below close the loop.
        if [ ! -s "\$READS" ]; then
            :
        elif [ "${o.match_aligner}" = "BLAST" ]; then
            makeblastdb -in "\$REFS" -dbtype nucl -out mdb > /dev/null 2>&1
            blastn -task ${o.blast_task} -query "\$READS" -db mdb -num_threads ${task.cpus} \\
                -dust ${o.blast_dust} -evalue ${o.blast_evalue} -max_target_seqs 7 \\
                -outfmt "6 qseqid sseqid pident length qstart qend sstart send sstrand bitscore" > hits.tsv 2>/dev/null
            blast_match.py --blast hits.tsv --query "\$READS" --out R\${round} \\
                --bin-map "\$REFS" \\
                --min-len ${o.blast_min_hsp_len} --min-pid ${o.blast_min_hsp_pid} ${incl}
            rm -f mdb.*
        else
            rammap -x sr --secondary yes -N 7 -c --filter-chimera \\
                --chimera-min-len ${o.blast_min_hsp_len} --chimera-min-pid ${o.blast_min_hsp_pid} \\
                -t ${task.cpus} "\$REFS" "\$READS" > hits.paf 2>/dev/null
            rammap_partition.py --paf hits.paf --query "\$READS" --out R\${round} \\
                --bin-map "\$REFS" \\
                --min-len ${o.blast_min_hsp_len} --min-pid ${o.blast_min_hsp_pid} ${incl} ${defer}
        fi
        for ext in match class nomatch chim; do [ -e R\${round}.\$ext ] || : > R\${round}.\$ext; done

        if [ ! -s R\${round}.class ]; then
            : > R\${round}.sort_stats.txt
            DONE=1; REASON=nothing_recruited
        else
            # ── SORT ────────────────────────────────────────────────────────────────────────────────
            RL=\$(awk '/^>/{next}{s+=length(\$0);n++} END{print (n>0?int(s/n):150)}' R\${round}.match)
            parse_sort_results.py R\${round}.class R\${round}.match R\${round} -C \$MRC -D \$MRP ${sg} ${bg} ${sec} ${collapse}
            [ -e R\${round}.txt ] && mv R\${round}.txt R\${round}.sort_stats.txt || true
            [ -e R\${round}.sort_stats.txt ] || : > R\${round}.sort_stats.txt
            rm -f R\${round}-*.fa.2 2>/dev/null || true

            # ── per gene: ALIGN -> REFINE ───────────────────────────────────────────────────────────
            for fa in R\${round}-*.fa; do
                [ -e "\$fa" ] || continue
                gene=\$(basename "\$fa" .fa); gene=\${gene#R\${round}-}
                # This gene's own reference from the prior round when it has one; in round 1 there is no
                # prior, so select this gene's record OUT of the panel rather than passing the panel.
                # Passing the whole panel meant sam_align_stats.py took ref_len from whichever record
                # sorted first — so reordering a panel of byte-identical sequences moved pct_called by
                # 12 points — and let a bin's reads align to a neighbour's coordinates (~0.2%).
                gref="\$REFS"
                prev=\$(( round - 1 ))
                if [ -e "R\${prev}-\${gene}.ref" ]; then
                    gref="R\${prev}-\${gene}.ref"
                else
                    # --subtype-counts: choose WHICH record in the bin on this sample's own read
                    # evidence, not on file order. rammap_partition.py wrote the tally in this same
                    # directory. Without it the bin always aligned against its first record, which
                    # build_panel emits RefSeq-first — canonical, but not necessarily the closest
                    # variant to this sample. On the arbovirus cohort that cost 77% of multi-record
                    # bins their best reference, and one JCV M segment went 41.5% -> 100% breadth
                    # once the right sibling was picked.
                    select_reference.py --fasta "\$REFS" --name "\${gene}" \\
                        --subtype-counts "R\${round}.subtypes.tsv" \\
                        -o "R\${round}-\${gene}.gref"
                    gref="R\${round}-\${gene}.gref"
                fi
                if [ "${o.gather_aligner}" = "BLAST" ]; then
                    makeblastdb -in "\$gref" -dbtype nucl -out adb > /dev/null 2>&1
                    blastn -task ${o.blast_task} -query "\$fa" -db adb -num_threads ${task.cpus} \\
                        -dust ${o.blast_dust} -evalue ${o.blast_evalue} -max_target_seqs 1 -max_hsps 1 \\
                        -outfmt "6 sstart send sstrand qseq sseq" 2>/dev/null \\
                        | blast_align_stats.py --ref "\$gref" -o R\${round}-\${gene}.json ${elo}
                    rm -f adb.*
                else
                    rammap -a -x sr --secondary no --eqx -t ${task.cpus} "\$gref" "\$fa" 2>/dev/null \\
                        | samtools view -F 4 -h - \\
                        | sam_align_stats.py --ref "\$gref" -o R\${round}-\${gene}.json ${elo} --term-window ${tw}
                fi
                # CUMULATIVE JSON history: this round's first, then every prior round's for this gene.
                # Cumulative is required, not an optimisation — a reference built from a single round's
                # JSON collapses on a round that recruits few reads, and that propagates into assembly.
                priors=\$( { ls R*-\${gene}.json 2>/dev/null | grep -v "^R\${round}-" | sort; } || true )
                keep=""
                [ "${o.gather_del_type}" = "REF" ] && keep="--keep-deleted \$gref"
                combine_align_stats.py --name "\$gene" --count-alt ${o.min_ca} --count-freq ${o.min_fa} \\
                    ${elo} ${del_ambig} \$keep ${pb} ${mrd} R\${round}-\${gene}.json \$priors > R\${round}-\${gene}.ref
            done

            # ── CLOSE: combined next ref + convergence ───────────────────────────────────────────────
            cat R\${round}-*.ref > next_refs.fasta 2>/dev/null || true
            [ -s next_refs.fasta ] || cp "\$REFS" next_refs.fasta

            # Two ways a round converges: the reference stopped moving, or nothing is left to recruit.
            # The second matters on shallow runs where round 1 consumes every read.
            if cmp -s "\$REFS" next_refs.fasta; then
                DONE=1; REASON=reference_stable
            elif [ ! -s R\${round}.nomatch ]; then
                DONE=1; REASON=nothing_left
            fi
            mv next_refs.fasta refs_R\${round}.fasta
            REFS=refs_R\${round}.fasta
            READS=R\${round}.nomatch
        fi

        t1=\$(date +%s)
        printf '%s\\t%s\\t%s\\t%s\\n' "\$round" "\$(( t1 - t0 ))" "\$DONE" "\$REASON" >> round_timings.tsv
        [ "\$DONE" -eq 1 ] && break
    done

    # Terminal reporting files, named here rather than in a downstream process. This task is the only
    # place that knows the round order, so it is the only place that can concatenate the per-round stats
    # correctly: a downstream `R*.sort_stats.txt` glob sorts lexically and would put R10 before R2,
    # which bites at max_rounds >= 10 (flu/sensitive sets exactly that).
    cp R\${LAST_ROUND}.nomatch nomatch.fa
    cp R1.chim chimeric.fa
    printf 'Gene\\tRead Patterns\\tRead Count\\n' > sorted_read_stats.txt
    for r in \$(seq 1 \$LAST_ROUND); do
        [ -e "R\${r}.sort_stats.txt" ] && cat "R\${r}.sort_stats.txt" >> sorted_read_stats.txt
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rammap: \$(rammap --version 2>/dev/null | grep -oE '[0-9][^ ]*' | head -1 || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    : > R1-TARGET.ref
    : > R1-TARGET.fa
    : > nomatch.fa
    : > chimeric.fa
    printf 'Gene\\tRead Patterns\\tRead Count\\n' > sorted_read_stats.txt
    printf 'round\\tseconds\\tconverged\\treason\\n1\\t0\\t1\\treference_stable\\n' > round_timings.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rammap: \$(rammap --version 2>/dev/null | grep -oE '[0-9][^ ]*' | head -1 || echo unknown)
    END_VERSIONS
    """
}
