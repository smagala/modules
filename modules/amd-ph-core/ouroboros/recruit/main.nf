// OUROBOROS_RECRUIT — one whole recruit round in a SINGLE process, driven by a topic feedback loop.
//
// The process reads a round state, runs the round (MATCH -> SORT -> per-gene ALIGN+REFINE ->
// convergence close), and emits the NEXT state. The calling subworkflow feeds that back in, so the
// loop lives there rather than here.
//
// The round is ONE process rather than five deliberately. An earlier decomposition into discrete
// MATCH/SORT/ALIGN/REFINE processes was correct on paper and lost badly on the cloud: the round's own
// compute is seconds, so per-step task overhead dominated. See the notebook's ADR-0013 for the trace
// data and for the R1-versus-tail split that is planned next.
//
// The round state is a TUPLE, staged into the instate/ layout by Nextflow rather than assembled by a
// prior process:
//   refs.fasta   combined reference matched this round
//   carry.fa     reads still to match — the prior round's nomatch, or the initial query
//   gene_refs/   prior per-gene refs, used for reference selection
//   jsons/       CUMULATIVE per-gene alignment JSONs. Cumulative is not an optimisation: building a
//                reference from only the current round's JSON collapses it on a round that recruits
//                few reads, which then propagates into assembly.
//
// Round 1 passes [] for gene_refs and jsons, so those directories simply do not exist — every read of
// them in the script below is already guarded. gene_refs/ and jsons/ are emitted as DIRECTORIES rather
// than as globs, because a barren round legitimately writes nothing into them and a `/*` output glob
// would fail the task on a missing match.
//
// The next state is emitted BOTH as a named output and onto the 'gather' topic. The topic is what
// closes the feedback loop, and it has to live here: only a process can publish to a topic, there is
// no operator that pushes a channel onto one. So this module does know the topic's name — an honest
// coupling, not something the subworkflow can take over. What the subworkflow does own is the seeding
// and the termination gate, which is why there is no seeding process and no `until` in here.
//
// Per-round artifacts that the assembly step collects are NORMAL (non-topic) outputs, and optional, so
// a barren round simply emits fewer of them rather than failing on a missing output.

process OUROBOROS_RECRUIT {
    tag "${meta.id}:R${round}"
    label 'process_high'

    container 'oamd-bio-ouroboros:1.0.0'

    input:
    tuple val(meta), path(refs,      stageAs: 'instate/refs.fasta'),
                     path(carry,     stageAs: 'instate/carry.fa'),
                     path(gene_refs, stageAs: 'instate/gene_refs'),
                     path(jsons,     stageAs: 'instate/jsons'),
                     val(round), val(done)   // `done` is the until gate; unused here
    val recruit_opts

    output:
    tuple val(meta), path("outstate/refs.fasta"), path("outstate/carry.fa"),
                     path("outstate/gene_refs"), path("outstate/jsons"),
                     env('RN'), env('DONE'), emit: state, topic: 'gather'
    tuple val(meta), path("R${round}-*.ref")           , val(round), emit: gene_refs  , optional: true
    tuple val(meta), path("R${round}-*.fa")            , val(round), emit: sorted     , optional: true
    tuple val(meta), path("R${round}.nomatch")         , val(round), emit: raw_nomatch
    tuple val(meta), path("R${round}.chim")            , val(round), emit: chimeric
    tuple val(meta), path("R${round}.sort_stats.txt")  , val(round), emit: stats
    path "versions.yml"                                , emit: versions

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
    def missing = required.findAll { !recruit_opts.containsKey(it) }
    if (missing) {
        error "OUROBOROS_RECRUIT: recruit_opts is missing required key(s): ${missing.join(', ')}"
    }

    def o     = recruit_opts
    def incl  = o.blast_incl_chim ? '--incl-chim' : ''
    // Round 1 applies the real recruitment floors; later rounds accept any read that matches, because
    // by then the reference has moved and a single read is meaningful evidence.
    def mrc   = (round as int) == 1 ? o.min_recruit_count : 1
    def mrp   = (round as int) == 1 ? o.min_recruit_patterns : 1
    def sg    = o.sort_groups ? "-P '${o.sort_groups}'" : ''
    def bg    = o.ban_groups  ? "-B '${o.ban_groups}'"  : ''
    // Secondary-target keep (co-infection / contamination): prefer a COVERAGE floor on the secondary's
    // own reference, which is invariant to how deep the primary is, over the fraction-of-primary gate,
    // which moves with the primary's depth. \$REFS and \$RL resolve in the script body.
    def sec   = (o.secondary_min_depth && (o.secondary_min_depth as float) > 0)
                    ? "--ref-fasta \$REFS --read-len \$RL --secondary-min-depth ${o.secondary_min_depth}"
                    : ((o.secondary_min_frac && (o.secondary_min_frac as float) > 0)
                        ? "--secondary-frac ${o.secondary_min_frac} --secondary-floor ${o.min_ca}"
                        : '')
    def elo       = o.skip_elongation ? '--skip-elongation' : ''
    def tw        = o.gather_position_bias ? 3 : 0
    def del_ambig = o.gather_del_type == 'NNN' ? '--delete-by-ambiguity' : ''
    def pb        = o.gather_position_bias ? '--position-bias-frac 0.9 --position-bias-min-span 5' : ''
    def mrd       = (o.gather_ref_floor && o.gather_del_type == 'REF' && o.min_column_coverage > 0)
                        ? "--min-ref-depth ${o.min_column_coverage}" : ''
    def rprev     = (round as int) - 1
    """
    set -e
    mkdir -p outstate/gene_refs outstate/jsons
    REFS=instate/refs.fasta
    READS=instate/carry.fa

    # ── MATCH ───────────────────────────────────────────────────────────────────────────────────────
    # A round can inherit no carry reads (a shallow run where round 1 recruited everything). rammap
    # spins indefinitely on a zero-byte query rather than returning empty, so the aligner must never
    # see one: skip MATCH and let the empty-class branch below close the loop.
    if [ ! -s \$READS ]; then
        :
    elif [ "${o.match_aligner}" = "BLAST" ]; then
        makeblastdb -in \$REFS -dbtype nucl -out mdb > /dev/null 2>&1
        blastn -task ${o.blast_task} -query \$READS -db mdb -num_threads ${task.cpus} \\
            -dust ${o.blast_dust} -evalue ${o.blast_evalue} -max_target_seqs 7 \\
            -outfmt "6 qseqid sseqid pident length qstart qend sstart send sstrand bitscore" > hits.tsv 2>/dev/null
        blast_match.py --blast hits.tsv --query \$READS --out R${round} \\
            --min-len ${o.blast_min_hsp_len} --min-pid ${o.blast_min_hsp_pid} ${incl}
        rm -f mdb.*
    else
        rammap -x sr --secondary yes -N 7 -c --filter-chimera \\
            --chimera-min-len ${o.blast_min_hsp_len} --chimera-min-pid ${o.blast_min_hsp_pid} \\
            -t ${task.cpus} \$REFS \$READS > hits.paf 2>/dev/null
        rammap_partition.py --paf hits.paf --query \$READS --out R${round} \\
            --min-len ${o.blast_min_hsp_len} --min-pid ${o.blast_min_hsp_pid} ${incl}
    fi
    for ext in match class nomatch chim; do [ -e R${round}.\$ext ] || : > R${round}.\$ext; done

    if [ ! -s R${round}.class ]; then
        cp \$REFS outstate/refs.fasta
        : > outstate/carry.fa
        DONE=1
    else
        # ── SORT ────────────────────────────────────────────────────────────────────────────────────
        RL=\$(awk '/^>/{next}{s+=length(\$0);n++} END{print (n>0?int(s/n):150)}' R${round}.match)
        parse_sort_results.py R${round}.class R${round}.match R${round} -C ${mrc} -D ${mrp} ${sg} ${bg} ${sec}
        [ -e R${round}.txt ] && mv R${round}.txt R${round}.sort_stats.txt || true
        rm -f R${round}-*.fa.2 2>/dev/null || true

        # ── per gene: ALIGN -> REFINE ───────────────────────────────────────────────────────────────
        for fa in R${round}-*.fa; do
            [ -e "\$fa" ] || continue
            gene=\$(basename "\$fa" .fa); gene=\${gene#R${round}-}
            gref="\$REFS"; [ -e "instate/gene_refs/R${rprev}-\${gene}.ref" ] && gref="instate/gene_refs/R${rprev}-\${gene}.ref"
            if [ "${o.gather_aligner}" = "BLAST" ]; then
                makeblastdb -in "\$gref" -dbtype nucl -out adb > /dev/null 2>&1
                blastn -task ${o.blast_task} -query "\$fa" -db adb -num_threads ${task.cpus} \\
                    -dust ${o.blast_dust} -evalue ${o.blast_evalue} -max_target_seqs 1 -max_hsps 1 \\
                    -outfmt "6 sstart send sstrand qseq sseq" 2>/dev/null \\
                    | blast_align_stats.py --ref "\$gref" -o R${round}-\${gene}.json ${elo}
                rm -f adb.*
            else
                rammap -a -x sr --secondary no --eqx -t ${task.cpus} "\$gref" "\$fa" 2>/dev/null \\
                    | samtools view -F 4 -h - \\
                    | sam_align_stats.py --ref "\$gref" -o R${round}-\${gene}.json ${elo} --term-window ${tw}
            fi
            # Cumulative JSON history: this round's JSON first, then priors R1..R{round-1} from instate.
            priors=\$( { ls instate/jsons/R*-\${gene}.json 2>/dev/null | sort; } || true )
            keep=""
            [ "${o.gather_del_type}" = "REF" ] && keep="--keep-deleted \$gref"
            combine_align_stats.py --name "\$gene" --count-alt ${o.min_ca} --count-freq ${o.min_fa} \\
                ${elo} ${del_ambig} \$keep ${pb} ${mrd} R${round}-\${gene}.json \$priors > R${round}-\${gene}.ref
        done

        # ── CLOSE: combined next ref + reference-stability convergence ───────────────────────────────
        cat R${round}-*.ref > outstate/refs.fasta 2>/dev/null || true
        [ -s outstate/refs.fasta ] || cp \$REFS outstate/refs.fasta
        # Two ways a round converges: the reference stopped moving, or nothing is left to recruit. The
        # second matters on shallow runs where round 1 consumes every read — carrying an empty nomatch
        # forward would leave the next round with a zero-byte query.
        if cmp -s \$REFS outstate/refs.fasta; then : > outstate/carry.fa; DONE=1
        elif [ ! -s R${round}.nomatch ]; then : > outstate/carry.fa; DONE=1
        else cp R${round}.nomatch outstate/carry.fa; DONE=0; fi
    fi

    # Populate the next state dir: this round's per-gene refs, and the cumulative JSON set.
    cp R${round}-*.ref outstate/gene_refs/ 2>/dev/null || true
    find instate/jsons -name '*.json' -exec cp {} outstate/jsons/ \\; 2>/dev/null || true
    cp R${round}-*.json outstate/jsons/ 2>/dev/null || true

    [ -e R${round}.sort_stats.txt ] || : > R${round}.sort_stats.txt
    RN=\$(( ${round} + 1 ))

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rammap: \$(rammap --version 2>/dev/null | grep -oE '[0-9][^ ]*' | head -1 || echo unknown)
    END_VERSIONS
    """
}
