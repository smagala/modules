#!/usr/bin/env nextflow

//
// OUROBOROS_ASSEMBLE — polish and call every bin, one task per bin.
//
// The polish loop turns inside OUROBOROS_POLISH, so this subworkflow contains no loop mechanics: no
// `.recurse()`, no bundled per-sample state, no bounded groupTuple re-aggregation, and no
// `nextflow.preview.recursion`. Each row of the input channel is one bin, and bins are an ordinary
// scatter — they run concurrently and converge independently.
//
// What that removes is worth naming, because it was all workaround. `.recurse()` is value-channel
// only, so the previous design bundled every gene of a sample into ONE recursion state, fanned out per
// gene each round, and re-aggregated under a bounded groupKey. That forced a ROUND BARRIER: every bin
// waited on the slowest bin at every iteration. A skip-converged filter was then added so a converged
// bin would stop re-aligning against an identical reference — a workaround for a barrier that only
// existed because of the bundling. With one task per bin, none of it is needed.
//

include { OUROBOROS_POLISH } from '../../../modules/amd-ph-core/ouroboros/polish/main'
include { OUROBOROS_CALL   } from '../../../modules/amd-ph-core/ouroboros/call/main'
include { OUROBOROS_REPORT } from '../../../modules/amd-ph-core/ouroboros/report/main'

workflow OUROBOROS_ASSEMBLE {

    take:
    ch_per_bin   // channel: [ val(meta), val(gene), path(sorted_reads), path(ref), path(deflation) ]
    align_opts   // value:   Map — aligner, sw_match, sw_mismatch, sw_gap_open, sw_gap_extend, term_ext,
                 //          ins_fold_freq, del_fold_freq, mut_min_strand_frac
    call_opts    // value:   Map — min_variant_count, min_variant_freq, min_ins_freq, min_del_freq,
                 //          min_variant_qual, min_column_coverage, min_confidence, sig_level,
                 //          auto_freq, min_ambig_freq, seg_numbers, sort_groups, caller, refchg, sor,
                 //          amend_indels, refchg_min_strand_frac
    max_iter     // value:   safety cap on polish iterations; convergence normally stops the loop first

    main:
    ch_versions = channel.empty()

    // Validate the option maps here so a typo fails with its key name rather than surfacing as a null
    // interpolated into a task script. Maps rather than ~26 positional values: at that arity positional
    // threading is unreadable and mis-binds silently on any reorder.
    def required_align = ['aligner', 'sw_match', 'sw_mismatch', 'sw_gap_open', 'sw_gap_extend',
                          'term_ext', 'ins_fold_freq', 'del_fold_freq', 'mut_min_strand_frac']
    def required_call  = ['min_variant_count', 'min_variant_freq', 'min_ins_freq', 'min_del_freq',
                          'min_variant_qual', 'min_column_coverage', 'min_confidence', 'sig_level',
                          'auto_freq', 'min_ambig_freq', 'seg_numbers', 'sort_groups', 'caller',
                          'refchg', 'sor', 'amend_indels', 'refchg_min_strand_frac']
    def missing_align = required_align.findAll { k -> !align_opts.containsKey(k) }
    def missing_call  = required_call.findAll  { k -> !call_opts.containsKey(k) }
    if (missing_align) {
        error "OUROBOROS_ASSEMBLE: align_opts is missing required key(s): ${missing_align.join(', ')}"
    }
    if (missing_call) {
        error "OUROBOROS_ASSEMBLE: call_opts is missing required key(s): ${missing_call.join(', ')}"
    }

    // One task per bin. Each converges on its own schedule.
    OUROBOROS_POLISH( ch_per_bin, align_opts, max_iter )
    ch_versions = ch_versions.mix( OUROBOROS_POLISH.out.versions.first() )

    // Route bins that aligned nothing away from calling. A 68-read bin whose reads the aligner cannot
    // anchor produces a zero-byte SAM, and variant calling on that fails — but the bin is a real
    // observation about the sample, so it is emitted on its own channel rather than silently dropped.
    ch_polished = OUROBOROS_POLISH.out.aligned.branch { _meta, _gene, _ref, _sam, reason ->
        callable    : reason != 'no_alignments'
        unalignable : true
    }

    OUROBOROS_CALL(
        ch_polished.callable.map { meta, gene, ref, sam, _reason -> [ meta, gene, ref, sam ] },
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
    consensus     = OUROBOROS_CALL.out.consensus       // channel: [ val(meta), val(gene), path(fasta) ]
    variants      = OUROBOROS_CALL.out.variants        // channel: [ val(meta), val(gene), path(txt) ]
    alleles       = OUROBOROS_CALL.out.alleles         // channel: [ val(meta), val(gene), path(txt) ]
    coverage      = OUROBOROS_CALL.out.coverage        // channel: [ val(meta), val(gene), path(txt) ]
    bam           = OUROBOROS_CALL.out.bam             // channel: [ val(meta), val(gene), path(bam) ]
    vcf           = OUROBOROS_CALL.out.vcf             // channel: [ val(meta), val(gene), path(vcf) ]
    amended       = OUROBOROS_CALL.out.amended         // channel: [ val(meta), val(gene), path(fa) ]
    amended_len   = OUROBOROS_CALL.out.amended_len     // channel: [ val(meta), val(gene), val(length) ]
    indels_folded = OUROBOROS_POLISH.out.indels_folded // channel: [ val(meta), val(gene), val(0|1) ]
    timings       = OUROBOROS_POLISH.out.timings       // channel: [ val(meta), val(gene), path(iter_timings.tsv) ]
    unalignable   = ch_polished.unalignable.map { meta, gene, _ref, _sam, reason -> [ meta, gene, reason ] } // channel: [ val(meta), val(gene), val(reason) ]
    versions      = ch_versions                        // channel: [ path(versions.yml) ]
}
