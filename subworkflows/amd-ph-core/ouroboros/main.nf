#!/usr/bin/env nextflow

//
// OUROBOROS — the complete recruit → refine → polish → call core, for a set of samples.
//
// Three stages, and the shape of each is set by what it compares:
//
//   OUROBOROS_GATHER          per SAMPLE — recruits reads into bins and refines a reference per bin
//   OUROBOROS_ASSEMBLE_STAGE  per BIN    — polishes and calls each bin independently
//   OUROBOROS_GENOMEREPORT    per SAMPLE — compares every genome of a sample against every other
//
// The middle stage is a plain scatter: OUROBOROS_GATHER emits one row per bin, and each becomes its own
// task that converges on its own schedule. The genomes screen has to come back to per-sample because it
// is inherently cross-bin — a contaminant is only identifiable relative to the sample's other genomes —
// so it waits on all of a sample's bins via groupTuple.
//
// Both loops turn inside their tasks. Nothing here iterates, and no state is written to disk between
// stages: the gather → assemble boundary travels as a channel.
//

include { OUROBOROS_GATHER         } from '../ouroboros_gather/main'
include { OUROBOROS_ASSEMBLE_STAGE } from '../ouroboros_assemble_stage/main'
include { OUROBOROS_GENOMEREPORT   } from '../../../modules/amd-ph-core/ouroboros/genomereport/main'

workflow OUROBOROS {

    take:
    ch_reads_refs       // channel: [ val(meta), path(reads_fasta), path(ref_set) ]
    ch_deflation        // channel: [ val(meta), path(deflation_index) ]
    recruit_opts        // value:   Map — see modules/amd-ph-core/ouroboros/recruit
    align_opts          // value:   Map — see modules/amd-ph-core/ouroboros/polish
    call_opts           // value:   Map — see modules/amd-ph-core/ouroboros/call
    max_rounds          // value:   safety cap on recruit rounds
    max_iter            // value:   safety cap on polish iterations
    amplicon_mode       // value:   boolean, force the reads-on-reference path for every bin
    stitch              // value:   boolean, project consensuses into the panel reference frame
    amplicon_min_depth  // value:   depth at or above which a reference position is called, not masked

    main:
    ch_versions = channel.empty()

    // ── recruit + refine, one task per sample ──
    OUROBOROS_GATHER( ch_reads_refs, ch_deflation, recruit_opts, max_rounds )
    ch_versions = ch_versions.mix( OUROBOROS_GATHER.out.versions )

    // ── polish + call, one task per BIN ──
    OUROBOROS_ASSEMBLE_STAGE(
        OUROBOROS_GATHER.out.per_gene,
        OUROBOROS_GATHER.out.ref_set,
        ch_deflation,
        align_opts,
        call_opts,
        max_iter,
        amplicon_mode,
        stitch,
        amplicon_min_depth
    )
    ch_versions = ch_versions.mix( OUROBOROS_ASSEMBLE_STAGE.out.versions )

    // ── genomes / contamination screen, back to one task per sample ──
    // Cross-bin by nature, so it gathers every bin's ref-pinned consensus for a sample before running.
    // groupTuple is safe here: these channels close normally, unlike inside a loop.
    ch_genome_report_in = OUROBOROS_ASSEMBLE_STAGE.out.refpinned
        .map { meta, _gene, fa -> [ meta, fa ] }
        .groupTuple()
        .join( OUROBOROS_GATHER.out.ref_set )
        .join( OUROBOROS_GATHER.out.sort_stats )

    OUROBOROS_GENOMEREPORT( ch_genome_report_in )
    ch_versions = ch_versions.mix( OUROBOROS_GENOMEREPORT.out.versions.first() )

    emit:
    genomes            = OUROBOROS_GENOMEREPORT.out.report          // channel: [ val(meta), path(genomes.tsv) ]
    consensus          = OUROBOROS_ASSEMBLE_STAGE.out.consensus     // channel: [ val(meta), val(gene), path(fasta) ]
    refpinned          = OUROBOROS_ASSEMBLE_STAGE.out.refpinned     // channel: [ val(meta), val(gene), path(fa) ]
    variants           = OUROBOROS_ASSEMBLE_STAGE.out.variants      // channel: [ val(meta), val(gene), path(txt) ]
    alleles            = OUROBOROS_ASSEMBLE_STAGE.out.alleles       // channel: [ val(meta), val(gene), path(txt) ]
    coverage           = OUROBOROS_ASSEMBLE_STAGE.out.coverage      // channel: [ val(meta), val(gene), path(txt) ]
    bam                = OUROBOROS_ASSEMBLE_STAGE.out.bam           // channel: [ val(meta), val(gene), path(bam) ]
    vcf                = OUROBOROS_ASSEMBLE_STAGE.out.vcf           // channel: [ val(meta), val(gene), path(vcf) ]
    amended            = OUROBOROS_ASSEMBLE_STAGE.out.amended       // channel: [ val(meta), val(gene), path(fa) ]
    stitched           = OUROBOROS_ASSEMBLE_STAGE.out.stitched      // channel: [ val(meta), val(gene), path(fa) ]
    relabelled         = OUROBOROS_ASSEMBLE_STAGE.out.relabelled    // channel: [ val(meta), val(gene), path(fa) ] — publish as <gene>.stitched.fa
    amplicon_bam       = OUROBOROS_ASSEMBLE_STAGE.out.amplicon_bam  // channel: [ val(meta), val(gene), path(bam), path(bai) ]
    amplicon_consensus = OUROBOROS_ASSEMBLE_STAGE.out.amplicon_consensus // channel: [ val(meta), val(gene), path(fa) ]
    amplicon_bed       = OUROBOROS_ASSEMBLE_STAGE.out.amplicon_bed  // channel: [ val(meta), val(gene), path(bed) ]
    chimeric           = OUROBOROS_GATHER.out.chimeric              // channel: [ val(meta), path(chimeric.fa) ]
    nomatch            = OUROBOROS_GATHER.out.nomatch               // channel: [ val(meta), path(nomatch.fa) ]
    sort_stats         = OUROBOROS_GATHER.out.sort_stats            // channel: [ val(meta), path(sorted_read_stats.txt) ]
    gather_timings     = OUROBOROS_GATHER.out.timings               // channel: [ val(meta), path(round_timings.tsv) ]
    collapse           = OUROBOROS_GATHER.out.collapse              // channel: [ val(meta), path(R*.collapse.tsv) ] — only when collapse_groups is set
    polish_timings     = OUROBOROS_ASSEMBLE_STAGE.out.timings       // channel: [ val(meta), val(gene), path(iter_timings.tsv) ]
    unalignable        = OUROBOROS_ASSEMBLE_STAGE.out.unalignable   // channel: [ val(meta), val(gene), val(reason) ] — bins the aligner could not anchor
    versions           = ch_versions                                // channel: [ path(versions.yml) ]
}
