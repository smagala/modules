#!/usr/bin/env nextflow

//
// OUROBOROS_ASSEMBLE_STAGE — the polish-and-call half of the pipeline, as a self-contained stage.
//
// Consumes the gather -> assemble boundary (per gene: sorted reads, refined reference, deflation index)
// plus the panel reference set, and runs OUROBOROS_ASSEMBLE followed by the two-way stitch gate.
// This is the per-bin unit: one invocation is one bin's independent recursion.
//
// The stitch gate exists so the reference-anchored read mapping runs ONLY where it can change the
// output. A gene needs the full reads-on-reference path when amplicon mode is on or when assembly
// folded an indel; otherwise the amended consensus already IS the stitch result and is relabelled
// rather than recomputed. That is the difference between mapping every gene's reads a second time and
// mapping almost none of them.
//

include { OUROBOROS_ASSEMBLE     } from '../ouroboros_assemble/main'
include { OUROBOROS_AMPLICONCALL } from '../../../modules/amd-ph-core/ouroboros/ampliconcall/main'
include { OUROBOROS_STITCH       } from '../../../modules/amd-ph-core/ouroboros/stitch/main'

workflow OUROBOROS_ASSEMBLE_STAGE {

    take:
    ch_per_gene         // channel: [ val(meta), val(gene), path(sorted_reads), path(refined_ref), path(deflation) ]
    ch_ref_set          // channel: [ val(meta), path(panel_ref_set) ]
    ch_deflation        // channel: [ val(meta), path(deflation_index) ]
    align_opts          // value:   Map, see OUROBOROS_ASSEMBLE
    call_opts           // value:   Map, see OUROBOROS_ASSEMBLE
    max_iter            // value:   safety cap on polishing rounds
    amplicon_mode       // value:   boolean, force the reads-on-reference path for every gene
    stitch              // value:   boolean, project the amended consensus into the reference frame
    amplicon_min_depth  // value:   depth at or above which a reference position is called rather than masked

    main:
    ch_versions = channel.empty()

    OUROBOROS_ASSEMBLE( ch_per_gene, align_opts, call_opts, max_iter )
    ch_versions = ch_versions.mix( OUROBOROS_ASSEMBLE.out.versions.first() )

    ch_sorted_by_gene = ch_per_gene.map { meta, gene, reads, _ref, _defl -> [ meta, gene, reads ] }

    ch_amplicon_bam       = channel.empty()
    ch_amplicon_consensus = channel.empty()
    ch_amplicon_bed       = channel.empty()
    ch_stitched           = channel.empty()
    ch_relabelled         = channel.empty()
    ch_refpinned          = channel.empty()

    if (amplicon_mode || stitch) {
        // Route each gene: run the full path only where it can change the output.
        ch_gene_route = OUROBOROS_ASSEMBLE.out.indels_folded
            .map { meta, gene, folded -> [ meta, gene, (amplicon_mode || (folded as int) != 0) ] }

        ch_run_genes      = ch_gene_route.filter { _meta, _gene, run ->  run }.map { meta, gene, _run -> [ meta, gene ] }
        ch_relabel_genes  = ch_gene_route.filter { _meta, _gene, run -> !run }.map { meta, gene, _run -> [ meta, gene ] }

        ch_amplicon_in = ch_sorted_by_gene
            .combine( ch_deflation, by: 0 )
            .combine( ch_ref_set, by: 0 )
            .map { meta, gene, reads, defl, refs -> [ meta, gene, reads, defl, refs ] }
            .join( ch_run_genes, by: [0, 1] )

        OUROBOROS_AMPLICONCALL( ch_amplicon_in, amplicon_min_depth )
        ch_amplicon_bam       = OUROBOROS_AMPLICONCALL.out.bam
        ch_amplicon_consensus = OUROBOROS_AMPLICONCALL.out.consensus
        ch_amplicon_bed       = OUROBOROS_AMPLICONCALL.out.bed
        ch_versions           = ch_versions.mix( OUROBOROS_AMPLICONCALL.out.versions.first() )

        if (stitch) {
            ch_stitch_in = OUROBOROS_ASSEMBLE.out.amended
                .join( OUROBOROS_AMPLICONCALL.out.bam, by: [0, 1] )
                .join( OUROBOROS_AMPLICONCALL.out.ref, by: [0, 1] )

            OUROBOROS_STITCH( ch_stitch_in, amplicon_min_depth )
            ch_stitched = OUROBOROS_STITCH.out.stitched
            ch_versions = ch_versions.mix( OUROBOROS_STITCH.out.versions.first() )

            // A full-length clean gene's stitched output is byte-identical to its amended consensus, so
            // it is relabelled rather than recomputed. Emitted as its own channel for the caller to
            // publish — a subworkflow should not reach for params.outdir to publish on its own.
            ch_relabelled = OUROBOROS_ASSEMBLE.out.amended.join( ch_relabel_genes, by: [0, 1] )

            // Reference-pinned consensus per gene = stitched (routed genes) + amended (relabelled genes).
            ch_refpinned = ch_stitched.mix( ch_relabelled )
        }
    }

    // When stitch is off there is no reference frame to pin to, so the reported consensus stands in.
    ch_refpinned_out = stitch ? ch_refpinned : OUROBOROS_ASSEMBLE.out.consensus

    emit:
    refpinned          = ch_refpinned_out                    // channel: [ val(meta), val(gene), path(fa) ]
    relabelled         = ch_relabelled                       // channel: [ val(meta), val(gene), path(fa) ] — publish as <gene>.stitched.fa
    consensus          = OUROBOROS_ASSEMBLE.out.consensus    // channel: [ val(meta), val(gene), path(fasta) ]
    variants           = OUROBOROS_ASSEMBLE.out.variants     // channel: [ val(meta), val(gene), path(txt) ]
    alleles            = OUROBOROS_ASSEMBLE.out.alleles      // channel: [ val(meta), val(gene), path(txt) ]
    coverage           = OUROBOROS_ASSEMBLE.out.coverage     // channel: [ val(meta), val(gene), path(txt) ]
    bam                = OUROBOROS_ASSEMBLE.out.bam          // channel: [ val(meta), val(gene), path(bam) ]
    vcf                = OUROBOROS_ASSEMBLE.out.vcf          // channel: [ val(meta), val(gene), path(vcf) ]
    amended            = OUROBOROS_ASSEMBLE.out.amended      // channel: [ val(meta), val(gene), path(fa) ]
    amplicon_bam       = ch_amplicon_bam                     // channel: [ val(meta), val(gene), path(bam), path(bai) ]
    amplicon_consensus = ch_amplicon_consensus               // channel: [ val(meta), val(gene), path(fa) ]
    amplicon_bed       = ch_amplicon_bed                     // channel: [ val(meta), val(gene), path(bed) ]
    stitched           = ch_stitched                         // channel: [ val(meta), val(gene), path(fa) ]
    timings            = OUROBOROS_ASSEMBLE.out.timings      // channel: [ val(meta), val(gene), path(iter_timings.tsv) ]
    unalignable        = OUROBOROS_ASSEMBLE.out.unalignable  // channel: [ val(meta), val(gene), val(reason) ]
    versions           = ch_versions                         // channel: [ path(versions.yml) ]
}
