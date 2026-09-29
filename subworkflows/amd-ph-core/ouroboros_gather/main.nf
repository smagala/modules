#!/usr/bin/env nextflow

//
// OUROBOROS_GATHER — the recruit/refine half of the pipeline, and the boundary it hands to assembly.
//
// The recruit loop turns INSIDE OUROBOROS_RECRUIT, one task per sample, so this subworkflow contains no
// loop mechanics at all: no topic, no `until`, no `.recurse()`. Samples are an ordinary scatter and
// converge independently.
//
// That is a correctness requirement, not a simplification. `until` closes the channel rather than
// dropping the matching item, so a topic feedback loop truncates every other sample the moment the
// first one converges — silently, because truncated recruitment merely yields a short reference.
// `filter` deadlocks instead, and `.recurse()` is value-channel only. See the notebook finding
// "until truncates per-bin scatter" and ADR-0013.
//
// What remains here is channel shaping: collect the loop's per-round artifacts into the gather ->
// assemble BOUNDARY — per gene, the sorted reads across all rounds, the HIGHEST-round refined
// reference, and the deflation index. That is the contract one assemble unit consumes, which is what
// lets assembly scatter per bin.
//
// (pirma split this across gather_recurse and gather_stage. The latter called no modules of its own —
// it only reshaped channels — so the two are merged rather than preserved as an indirection.)
//

include { OUROBOROS_RECRUIT } from '../../../modules/amd-ph-core/ouroboros/recruit/main'

workflow OUROBOROS_GATHER {

    take:
    ch_reads_refs   // channel: [ val(meta), path(reads_fasta), path(ref_set) ]
    ch_deflation    // channel: [ val(meta), path(deflation_index) ]
    recruit_opts    // value:   Map of recruit settings; see modules/amd-ph-core/ouroboros/recruit
    max_rounds      // value:   safety cap on recruit rounds

    main:
    ch_versions = channel.empty()

    if (recruit_opts.gather_aligner == 'NONE') {
        error "OUROBOROS_GATHER: gather_aligner = 'NONE' produces no refined references, so assembly " +
              "would have nothing to polish."
    }

    OUROBOROS_RECRUIT( ch_reads_refs, recruit_opts, max_rounds )
    ch_versions = ch_versions.mix( OUROBOROS_RECRUIT.out.versions.first() )

    // ── the gather -> assemble boundary ──
    // Per-gene sorted reads across every round; the round prefix is stripped so the key is the gene.
    ch_sorted_by_gene = OUROBOROS_RECRUIT.out.sorted
        .transpose()
        .map { meta, fa -> [ meta, fa.baseName.replaceAll(/^R\d+-/, ''), fa ] }
        .groupTuple( by: [0, 1] )

    // Each gene's HIGHEST-round reference: the one every prior round fed into.
    ch_refs_by_gene = OUROBOROS_RECRUIT.out.gene_refs
        .transpose()
        .map { meta, ref ->
            def gene  = ref.baseName.replaceAll(/^R\d+-/, '')
            def round = (ref.baseName =~ /^R(\d+)-/)[0][1] as int
            [ meta, gene, round, ref ]
        }
        .groupTuple( by: [0, 1] )
        .map { meta, gene, rounds, refs -> [ meta, gene, refs[ rounds.indexOf( rounds.max() ) ] ] }

    // One row per BIN. This is the scatter unit for assembly.
    ch_per_gene = ch_sorted_by_gene
        .join( ch_refs_by_gene, by: [0, 1] )
        .combine( ch_deflation, by: 0 )
        .map { meta, gene, fas, ref, deflation -> [ meta, gene, fas, ref, deflation ] }

    emit:
    per_gene = ch_per_gene                                                // channel: [ val(meta), val(gene), path(sorted_fas), path(refined_ref), path(deflation) ]
    ref_set  = ch_reads_refs.map { meta, _reads, refs -> [ meta, refs ] }  // channel: [ val(meta), path(panel_ref_set) ]
    chimeric   = OUROBOROS_RECRUIT.out.chimeric                           // channel: [ val(meta), path(chimeric.fa) ]
    nomatch    = OUROBOROS_RECRUIT.out.nomatch                            // channel: [ val(meta), path(nomatch.fa) ] — the unrecruited fraction behind every reported depth
    sort_stats = OUROBOROS_RECRUIT.out.sort_stats                         // channel: [ val(meta), path(sorted_read_stats.txt) ]
    timings  = OUROBOROS_RECRUIT.out.timings                              // channel: [ val(meta), path(round_timings.tsv) ]
    collapse = OUROBOROS_RECRUIT.out.collapse                             // channel: [ val(meta), path(R*.collapse.tsv) ] — only when collapse_groups is set
    versions = ch_versions                                              // channel: [ path(versions.yml) ]
}
