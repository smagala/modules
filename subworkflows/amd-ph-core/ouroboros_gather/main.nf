#!/usr/bin/env nextflow

//
// OUROBOROS_GATHER — the recruit/refine half of the pipeline, and the boundary it hands to assembly.
//
// The loop is a TOPIC FEEDBACK loop rather than `.recurse()`, and it is owned HERE rather than by the
// modules. The seed state is a plain channel map — round 1 is just the panel reference set, the query
// reads, and two empty directory slots — so there is no seeding *process*: Nextflow's input staging
// builds the instate/ layout, and a task that only ran `mkdir` and `cp` would be pure overhead.
// OUROBOROS_RECRUIT publishes its next state onto the 'gather' topic, and `until { round > max ||
// done }` drops the converged item to stop the loop.
//
// The topic name necessarily lives in the module, because only a process can publish to a topic — no
// operator pushes a channel onto one. What lives HERE is the seeding and the termination gate, so the
// module can be run for a single round, or tested, without a loop around it.
//
// Per-round artifacts come off OUROBOROS_RECRUIT's normal (non-topic) outputs and are collected after
// the loop closes.
//
// This subworkflow also shapes the gather -> assemble BOUNDARY: per gene, the sorted reads across all
// rounds plus the HIGHEST-round refined reference plus the deflation index. That is the contract one
// assemble run consumes, which is what allows assembly to be dispatched per bin.
//
// (pirma split this across gather_recurse and gather_stage. The latter called no modules of its own —
// it only reshaped channels — so the two are merged here rather than preserved as an indirection.)
//

include { OUROBOROS_RECRUIT       } from '../../../modules/amd-ph-core/ouroboros/recruit/main'
include { OUROBOROS_GATHERREPORTS } from '../../../modules/amd-ph-core/ouroboros/gatherreports/main'

workflow OUROBOROS_GATHER {

    take:
    ch_reads_refs   // channel: [ val(meta), path(reads_fasta), path(ref_set) ]
    ch_deflation    // channel: [ val(meta), path(deflation_index) ]
    recruit_opts    // value:   Map of recruit settings; see modules/amd-ph-core/ouroboros/recruit
    max_rounds      // value:   safety cap on recruit rounds

    main:
    ch_versions = Channel.empty()

    if (recruit_opts.gather_aligner == 'NONE') {
        error "OUROBOROS_GATHER: gather_aligner = 'NONE' produces no refined references, so assembly " +
              "would have nothing to polish."
    }

    // Round-1 state, as data rather than as a task. The two empty lists are gene_refs/ and jsons/:
    // there is no prior round to carry, and every read of those directories in the module is guarded.
    ch_seed = ch_reads_refs.map { meta, reads, refs -> [ meta, refs, reads, [], [], 1, 0 ] }

    // The loop: seed, then RECRUIT's own emissions fed back through the topic, truncated by `until`.
    ch_loop = ch_seed
        .mix( channel.topic('gather') )
        .until { meta, refs, carry, gene_refs, jsons, round, done ->
            (round as int) > (max_rounds as int) || (done as int) == 1
        }

    OUROBOROS_RECRUIT( ch_loop, recruit_opts )

    ch_versions = ch_versions.mix( OUROBOROS_RECRUIT.out.versions.first() )

    // ── collect per-round artifacts (these channels close when the loop stops) ──
    ch_sorted_reads = OUROBOROS_RECRUIT.out.sorted
        .flatMap { meta, fas, round -> (fas instanceof List ? fas : [fas]).collect { [ meta, it ] } }
        .groupTuple()

    ch_refined_refs = OUROBOROS_RECRUIT.out.gene_refs
        .flatMap { meta, refs, round -> (refs instanceof List ? refs : [refs]).collect { [ meta, it ] } }
        .groupTuple()

    // chimeric.fa is round 1's chimeras; nomatch.fa is the LAST round's raw nomatch; the sort stats are
    // every round's, concatenated in round order.
    ch_chimeric = OUROBOROS_RECRUIT.out.chimeric
        .filter { meta, chim, round -> (round as int) == 1 }
        .map    { meta, chim, round -> [ meta, chim ] }

    ch_nomatch = OUROBOROS_RECRUIT.out.raw_nomatch
        .map { meta, nomatch, round -> [ meta, [ (round as int), nomatch ] ] }
        .groupTuple()
        .map { meta, pairs -> [ meta, pairs.max { it[0] }[1] ] }

    ch_sort_stats = OUROBOROS_RECRUIT.out.stats
        .map { meta, stats, round -> [ meta, [ (round as int), stats ] ] }
        .groupTuple()
        .map { meta, pairs -> [ meta, pairs.sort { it[0] }.collect { it[1] } ] }

    OUROBOROS_GATHERREPORTS( ch_chimeric.join( ch_nomatch ).join( ch_sort_stats ) )

    // ── the gather -> assemble boundary ──
    ch_sorted_by_gene = ch_sorted_reads
        .transpose()
        .map { meta, fa -> [ meta, fa.baseName.replaceAll(/^R\d+-/, ''), fa ] }
        .groupTuple( by: [0, 1] )

    // Take each gene's HIGHEST-round reference: that is the one every prior round fed into.
    ch_refs_by_gene = ch_refined_refs
        .transpose()
        .map { meta, ref ->
            def gene  = ref.baseName.replaceAll(/^R\d+-/, '')
            def round = (ref.baseName =~ /^R(\d+)-/)[0][1] as int
            [ meta, gene, round, ref ]
        }
        .groupTuple( by: [0, 1] )
        .map { meta, gene, rounds, refs -> [ meta, gene, refs[ rounds.indexOf( rounds.max() ) ] ] }

    ch_per_gene = ch_sorted_by_gene
        .join( ch_refs_by_gene, by: [0, 1] )
        .combine( ch_deflation, by: 0 )
        .map { meta, gene, fas, ref, deflation -> [ meta, gene, fas, ref, deflation ] }

    emit:
    per_gene = ch_per_gene                                                   // channel: [ val(meta), val(gene), path(sorted_fas), path(refined_ref), path(deflation) ]
    ref_set  = ch_reads_refs.map { meta, reads, refs -> [ meta, refs ] }     // channel: [ val(meta), path(panel_ref_set) ]
    reports  = OUROBOROS_GATHERREPORTS.out.reports                           // channel: [ val(meta), path(chimeric), path(nomatch), path(sorted_read_stats) ]
    versions = ch_versions                                                   // channel: [ path(versions.yml) ]
}
