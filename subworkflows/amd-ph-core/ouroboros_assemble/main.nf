#!/usr/bin/env nextflow

//
// OUROBOROS_ASSEMBLE — iterative assembly polishing plus terminal variant calling.
//
// OUROBOROS_ASSEMBLE_ROUND is one polishing round expressed as native Nextflow recursion. The recurse
// state is ONE value per sample bundling ALL genes, because .recurse() takes a value channel; the round
// fans out per gene through OUROBOROS_ALIGN, then re-aggregates under a BOUNDED groupKey. The loop runs
// until every gene's reference reaches a fixpoint or the iteration cap is hit, and the best-scoring
// round's (ref, sam) pair feeds OUROBOROS_CALL.
//
// Recurse state, identical shape in and out: [ meta, genes, round, converged, align_opts ]
//   genes — one element per gene:
//     [ gene, align_ref, reads, done_ref, done_sam, converged, score, folded ]
//       align_ref           reference to align against THIS round (the prior round's amended ref)
//       reads               per-gene INFLATED reads, inflated once in the seed and carried unchanged,
//                           so ALIGNSPLIT only ever partitions and never re-inflates
//       done_ref, done_sam  the (ref, sam) actually aligned THIS round; what CALL consumes
//       converged           1 iff THIS gene reached a fixpoint
//       score               this round's alignment score, carried for the plateau test (long: a sum of
//                           per-read scores overflows int on very deep libraries)
//       folded              1 iff an indel was folded THIS round or ANY prior round (OR-accumulated)
//   converged (state)       1 iff ALL genes converged this round
//   align_opts              alignment settings, constant across rounds. These ride inside the state
//                           because .recurse() requires the take and emit shapes to match exactly, so a
//                           second `take:` is not available — and reaching for `params.*` inside the
//                           round would defeat the point of threading settings as inputs.
//
// Requires `nextflow.preview.recursion = true` in the entry script. A nextflow.config setting is
// ignored, and the flag is still required on 25.10.4 despite the docs describing it as deprecated.
//

include { OUROBOROS_INFLATE } from '../../../modules/amd-ph-core/ouroboros/inflate/main'
include { OUROBOROS_ALIGN   } from '../ouroboros_align/main'
include { OUROBOROS_CALL    } from '../../../modules/amd-ph-core/ouroboros/call/main'
include { OUROBOROS_REPORT  } from '../../../modules/amd-ph-core/ouroboros/report/main'

//
// One polishing round. Kept as a named workflow because .recurse() operates on a workflow.
//
workflow OUROBOROS_ASSEMBLE_ROUND {

    take:
    ch_state     // channel: [ val(meta), val(genes), val(round), val(converged), val(align_opts) ]

    main:
    // align_opts rides INSIDE the recursion state rather than arriving as a second `take:`.
    // .recurse() requires the take and emit shapes to match exactly, so a second take argument is not
    // available; carrying the settings as a constant element of the state keeps the contract intact and
    // keeps this subworkflow free of any `params.*` access.
    ch_opts = ch_state.map { meta, genes, round, converged, align_opts -> [ [ meta, round ], align_opts ] }
    // Fan out per gene, carrying the round.
    //
    // CRITICAL: under .recurse() several rounds' values are in flight at once, so EVERY per-gene join is
    // keyed by [meta, gene, round] and every re-aggregation by [meta, round]. Keying on [meta, gene]
    // alone pairs values across rounds by arrival order and silently corrupts the converged/score state.
    //
    // Skip-converged: only NOT-yet-converged genes are aligned each round; converged genes are carried
    // through frozen (ch_carried below). A converged gene would otherwise re-align against an identical
    // reference every round, so a thin contaminant would keep doing exact Smith-Waterman until the
    // DEEPEST gene converged. The terminal still picks each gene's best round from the rounds it ran.
    ch_gene = ch_state.flatMap { meta, genes, round, converged, align_opts ->
        genes.findAll { g -> ((g.size() > 5 ? g[5] : 0) as int) == 0 }
             .collect { g -> [ meta, g[0], round, g[2], g[1] ] }   // [ meta, gene, round, reads, align_ref ]
    }

    // The settings are constant across rounds, so first() gives OUROBOROS_ALIGN the value channels it
    // needs without re-deriving them per round.
    ch_align_opts = ch_state.map { meta, genes, round, converged, align_opts -> align_opts }.first()

    OUROBOROS_ALIGN(
        ch_gene,
        ch_align_opts.map { it.aligner },
        ch_align_opts.map { it.chunk_reads },
        ch_align_opts.map { it.sw_match },
        ch_align_opts.map { it.sw_mismatch },
        ch_align_opts.map { it.sw_gap_open },
        ch_align_opts.map { it.sw_gap_extend },
        ch_align_opts.map { it.term_ext },
        ch_align_opts.map { it.ins_fold_freq },
        ch_align_opts.map { it.del_fold_freq },
        ch_align_opts.map { it.mut_min_strand_frac }
    )

    // Re-join this round's result with the carried reads, prior score, and prior fold flag.
    ch_reads = ch_state.flatMap { meta, genes, round, converged, align_opts ->
        genes.collect { g -> [ [ meta, g[0], round ], g[2] ] }
    }
    ch_prev_score = ch_state.flatMap { meta, genes, round, converged, align_opts ->
        genes.collect { g -> [ [ meta, g[0], round ], ((g.size() > 6 ? g[6] : -1) as long) ] }
    }
    ch_prev_folded = ch_state.flatMap { meta, genes, round, converged, align_opts ->
        genes.collect { g -> [ [ meta, g[0], round ], ((g.size() > 7 ? g[7] : 0) as int) ] }
    }

    ch_result = OUROBOROS_ALIGN.out.sam
        .map { meta, gene, round, ref, sam, score -> [ [ meta, gene, round ], ref, sam, (score as long) ] }
        .join( OUROBOROS_ALIGN.out.refined.map { meta, gene, round, next_ref, conv, folded ->
            [ [ meta, gene, round ], next_ref, (conv as int), (folded as int) ] } )
        .join( ch_reads )
        .join( ch_prev_score )
        .join( ch_prev_folded )
        .map { key, ref, sam, score, next_ref, fixpoint, this_fold, reads, prev_score, prev_fold ->
            def round = key[2] as int
            // Convergence is the byte fixpoint OR a score plateau. The plateau matters for samples whose
            // reference jitters between two states without ever byte-fixpointing.
            def plateau = (round > 1 && score <= prev_score) ? 1 : 0
            def conv    = (fixpoint == 1 || plateau == 1) ? 1 : 0
            def folded  = (this_fold == 1 || prev_fold == 1) ? 1 : 0
            [ [ key[0], round ], [ key[1], next_ref, reads, ref, sam, conv, score, folded ] ]
        }

    // Converged genes carried forward untouched, so the bundle still holds all n_genes.
    ch_carried = ch_state.flatMap { meta, genes, round, converged, align_opts ->
        genes.findAll { g -> ((g.size() > 5 ? g[5] : 0) as int) == 1 }
             .collect { g -> [ [ meta, round ], g ] }
    }

    // Bounded re-aggregation per [meta, round]: collect exactly n_genes elements and rebuild the state.
    ch_n_genes = ch_state.map { meta, genes, round, converged, align_opts -> [ [ meta, round ], genes.size() ] }

    next_state = ch_result.mix( ch_carried )
        .combine( ch_n_genes, by: 0 )
        .map { key, element, n -> [ groupKey( key, n ), element ] }
        .groupTuple()
        .join( ch_opts )
        .map { key, elements, align_opts ->
            def all_converged = elements.every { ((it.size() > 5 ? it[5] : 0) as int) == 1 } ? 1 : 0
            [ key[0], elements, (key[1] as int) + 1, all_converged, align_opts ]
        }

    emit:
    next_state   // channel: [ val(meta), val(genes), val(round), val(converged), val(align_opts) ]
}

workflow OUROBOROS_ASSEMBLE {

    take:
    ch_in        // channel: [ val(meta), val(gene), path(sorted_reads), path(ref), path(deflation_index) ]
    align_opts   // value:   Map — aligner, chunk_reads, sw_match, sw_mismatch, sw_gap_open,
                 //          sw_gap_extend, term_ext, ins_fold_freq, del_fold_freq, mut_min_strand_frac
    call_opts    // value:   Map — min_variant_count, min_variant_freq, min_ins_freq, min_del_freq,
                 //          min_variant_qual, min_column_coverage, min_confidence, sig_level,
                 //          auto_freq, min_ambig_freq, seg_numbers, sort_groups, caller, refchg, sor,
                 //          amend_indels, refchg_min_strand_frac
    max_iter     // value:   safety cap on polishing rounds; convergence normally stops the loop first

    main:
    ch_versions = Channel.empty()

    // Settings arrive as Maps rather than ~25 positional values: positional threading at that arity is
    // unreadable and silently mis-binds on reorder. Validate the keys up front so a typo fails loudly
    // here instead of surfacing as a null inside a task script.
    def required_align = ['aligner', 'chunk_reads', 'sw_match', 'sw_mismatch', 'sw_gap_open',
                          'sw_gap_extend', 'term_ext', 'ins_fold_freq', 'del_fold_freq',
                          'mut_min_strand_frac']
    def required_call  = ['min_variant_count', 'min_variant_freq', 'min_ins_freq', 'min_del_freq',
                          'min_variant_qual', 'min_column_coverage', 'min_confidence', 'sig_level',
                          'auto_freq', 'min_ambig_freq', 'seg_numbers', 'sort_groups', 'caller',
                          'refchg', 'sor', 'amend_indels', 'refchg_min_strand_frac']
    def missing_align = required_align.findAll { !align_opts.containsKey(it) }
    def missing_call  = required_call.findAll  { !call_opts.containsKey(it) }
    if (missing_align) {
        error "OUROBOROS_ASSEMBLE: align_opts is missing required key(s): ${missing_align.join(', ')}"
    }
    if (missing_call) {
        error "OUROBOROS_ASSEMBLE: call_opts is missing required key(s): ${missing_call.join(', ')}"
    }

    // Inflate ONCE per gene. The reads are round-invariant — only the reference changes — so inflating
    // inside the recursion would repeat the same expansion every round.
    OUROBOROS_INFLATE( ch_in.map { meta, gene, reads, ref, defl -> [ meta, gene, reads, defl ] } )
    ch_versions = ch_versions.mix( OUROBOROS_INFLATE.out.versions.first() )

    ch_seed_in = ch_in.map { meta, gene, reads, ref, defl -> [ [ meta, gene ], ref ] }
        .join( OUROBOROS_INFLATE.out.fastq.map { meta, gene, fastq -> [ [ meta, gene ], fastq ] } )
        .map { key, ref, fastq -> [ key[0], key[1], fastq, ref ] }

    // Seed: bundle every gene of a sample into one state value. done_ref seeds to the gather reference;
    // done_sam is a placeholder until round 1 fills it.
    ch_seed = ch_seed_in
        .map { meta, gene, fastq, ref -> [ meta, [ gene, ref, fastq, ref, [], 0, -1, 0 ] ] }
        .groupTuple( by: 0 )
        .map { meta, genes -> [ meta, genes, 1, 0, align_opts ] }
        .first()

    OUROBOROS_ASSEMBLE_ROUND
        .recurse( ch_seed )
        .until { state -> (state[3] as int) == 1 || (state[2] as int) > (max_iter as int) }

    // Terminal: per gene, the HIGHEST-scoring round's (ref, sam) feeds CALL. recurse emits every round,
    // so picking the max-score round handles both a plateau and an outright score decrease without
    // threading a "previous best" through the state.
    ch_terminal = OUROBOROS_ASSEMBLE_ROUND.out
        .flatMap { meta, genes, round, converged, align_opts ->
            genes.collect { g -> [ [ meta, g[0] ], g[3], g[4], (g[6] as long), (g.size() > 7 ? (g[7] as int) : 0) ] }
        }
        .groupTuple()
        .map { key, done_refs, done_sams, scores, foldeds ->
            def i = scores.indexOf( scores.max() )
            // folded is OR-accumulated, so ANY round carrying 1 means an indel folded for this gene
            def folded = foldeds.any { it == 1 } ? 1 : 0
            [ key[0], key[1], done_refs[i], done_sams[i], folded ]
        }

    ch_indels_folded = ch_terminal.map { meta, gene, done_ref, done_sam, folded -> [ meta, gene, folded ] }

    OUROBOROS_CALL(
        ch_terminal.map { meta, gene, done_ref, done_sam, folded -> [ meta, gene, done_ref, done_sam ] },
        call_opts.min_variant_count,
        call_opts.min_variant_freq,
        call_opts.min_ins_freq,
        call_opts.min_del_freq,
        call_opts.min_variant_qual,
        call_opts.min_column_coverage,
        call_opts.min_confidence,
        call_opts.sig_level,
        call_opts.auto_freq,
        call_opts.min_ambig_freq,
        call_opts.seg_numbers,
        call_opts.sort_groups,
        call_opts.caller,
        call_opts.refchg,
        call_opts.sor,
        call_opts.amend_indels,
        call_opts.refchg_min_strand_frac
    )
    ch_versions = ch_versions.mix( OUROBOROS_CALL.out.versions.first() )

    // Diagnostic PDFs. pairing is paired-end only, so the remainder join fills it with [] when absent.
    ch_report_in = OUROBOROS_CALL.out.coverage
        .join( OUROBOROS_CALL.out.variants, by: [0, 1] )
        .join( OUROBOROS_CALL.out.alleles , by: [0, 1] )
        .join( OUROBOROS_CALL.out.pairing , by: [0, 1], remainder: true )
        .map { meta, gene, coverage, variants, alleles, pairing ->
            [ meta, gene, coverage, variants, alleles, pairing ?: [] ]
        }
    OUROBOROS_REPORT(
        ch_report_in,
        call_opts.min_variant_qual,
        call_opts.min_variant_freq,
        call_opts.min_column_coverage,
        call_opts.min_confidence
    )
    ch_versions = ch_versions.mix( OUROBOROS_REPORT.out.versions.first() )

    emit:
    consensus     = OUROBOROS_CALL.out.consensus   // channel: [ val(meta), val(gene), path(fasta) ]
    variants      = OUROBOROS_CALL.out.variants    // channel: [ val(meta), val(gene), path(txt) ]
    alleles       = OUROBOROS_CALL.out.alleles     // channel: [ val(meta), val(gene), path(txt) ]
    coverage      = OUROBOROS_CALL.out.coverage    // channel: [ val(meta), val(gene), path(txt) ]
    bam           = OUROBOROS_CALL.out.bam         // channel: [ val(meta), val(gene), path(bam) ]
    vcf           = OUROBOROS_CALL.out.vcf         // channel: [ val(meta), val(gene), path(vcf) ]
    amended       = OUROBOROS_CALL.out.amended     // channel: [ val(meta), val(gene), path(fa) ]
    amended_len   = OUROBOROS_CALL.out.amended_len // channel: [ val(meta), val(gene), val(length) ]
    indels_folded = ch_indels_folded               // channel: [ val(meta), val(gene), val(0|1) ]
    versions      = ch_versions                    // channel: [ path(versions.yml) ]
}
