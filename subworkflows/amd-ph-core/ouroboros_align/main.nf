#!/usr/bin/env nextflow

//
// OUROBOROS_ALIGN — one assembly-round alignment as a scatter-gather
// (ALIGNSPLIT -> ALIGNCHUNK -> ALIGNMERGE).
//
// ALIGNSPLIT partitions the per-gene, already-inflated reads into chunks (one ALIGNCHUNK task each, so
// the work packs across cores locally or scales across cloud nodes). A BOUNDED groupKey gather then
// feeds ALIGNMERGE, which produces the winnowed SAM, the assembly score, the pileup JSON, and the
// refined reference the next round consumes.
//
// Keying discipline (the gather_recurse lesson): under .recurse() multiple rounds flow through these
// channels concurrently, so EVERY per-chunk join is keyed by [meta, gene, round] and the gather is
// BOUNDED — an unbounded groupTuple hangs the recursion, because it waits for a channel close that
// never comes.
//

include { OUROBOROS_ALIGNSPLIT } from '../../../modules/amd-ph-core/ouroboros/alignsplit/main'
include { OUROBOROS_ALIGNCHUNK } from '../../../modules/amd-ph-core/ouroboros/alignchunk/main'
include { OUROBOROS_ALIGNMERGE } from '../../../modules/amd-ph-core/ouroboros/alignmerge/main'

workflow OUROBOROS_ALIGN {

    take:
    ch_gene              // channel: [ val(meta), val(gene), val(round), path(reads), path(align_ref) ]
    aligner              // value:   'RAMMAP' | 'LOCALSW'
    chunk_reads          // value:   reads per chunk above which the scatter fans out
    sw_match             // value:   Smith-Waterman match score
    sw_mismatch          // value:   Smith-Waterman mismatch penalty
    sw_gap_open          // value:   gap-open penalty, inclusive of the first extension (IRMA convention)
    sw_gap_extend        // value:   gap-extension penalty
    term_ext             // value:   reference-terminus end-bonus extension for local_sw
    ins_fold_freq        // value:   insertion frequency at or above which an indel folds into the ref
    del_fold_freq        // value:   deletion frequency at or above which an indel folds into the ref
    mut_min_strand_frac  // value:   min strand fraction for a reference-changing base (0 disables)

    main:
    ch_versions = Channel.empty()

    // ALIGNSPLIT emits no versions: it only partitions a FASTQ with coreutils.
    OUROBOROS_ALIGNSPLIT( ch_gene, aligner, chunk_reads )

    // chunk-file list -> one row per chunk, carrying the chunk id parsed from its filename
    ch_chunks = OUROBOROS_ALIGNSPLIT.out.chunks
        .transpose()
        .map { meta, gene, round, ref, chunk ->
            def chunk_id = chunk.getName().replaceAll(/^chunk\./, '').replaceAll(/\.fastq$/, '')
            [ meta, gene, round, chunk_id, ref, chunk ]
        }

    OUROBOROS_ALIGNCHUNK(
        ch_chunks,
        aligner,
        sw_match,
        sw_mismatch,
        sw_gap_open,
        sw_gap_extend,
        term_ext
    )
    ch_versions = ch_versions.mix( OUROBOROS_ALIGNCHUNK.out.versions.first() )

    // BOUNDED gather: groupKey(key, nchunk) emits as soon as N chunks arrive, so the recursion never
    // waits on a channel close.
    ch_nchunk = OUROBOROS_ALIGNSPLIT.out.nchunk
        .map { meta, gene, round, n -> [ [ meta, gene, round ], (n as int) ] }

    ch_align_ref = OUROBOROS_ALIGNSPLIT.out.chunks
        .map { meta, gene, round, ref, chunks -> [ [ meta, gene, round ], ref ] }

    ch_merge_in = OUROBOROS_ALIGNCHUNK.out.sam
        .map { meta, gene, round, sam -> [ [ meta, gene, round ], sam ] }
        .combine( ch_nchunk, by: 0 )
        .map { key, sam, n -> [ groupKey( key, n ), sam ] }
        .groupTuple()
        .join( ch_align_ref )
        .map { key, sams, ref -> [ key[0], key[1], key[2], ref, sams ] }

    OUROBOROS_ALIGNMERGE(
        ch_merge_in,
        aligner,
        ins_fold_freq,
        del_fold_freq,
        mut_min_strand_frac
    )
    ch_versions = ch_versions.mix( OUROBOROS_ALIGNMERGE.out.versions.first() )

    emit:
    sam      = OUROBOROS_ALIGNMERGE.out.sam      // channel: [ val(meta), val(gene), val(round), path(ref), path(sam), val(score) ]
    refined  = OUROBOROS_ALIGNMERGE.out.refined  // channel: [ val(meta), val(gene), val(round), path(next_ref), val(converged), val(indels_folded) ]
    pileup   = OUROBOROS_ALIGNMERGE.out.pileup   // channel: [ val(meta), val(gene), val(round), path(json) ]
    versions = ch_versions                       // channel: [ path(versions.yml) ]
}
