#!/usr/bin/env nextflow

//
// BBDUK_TRIM_2PASS: double primer trim on unaligned FASTQs
//

include { BBMAP_BBDUK as BBDUK_TRIM_1 } from "../../../modules/nf-core/bbmap/bbduk/main"
include { BBMAP_BBDUK as BBDUK_TRIM_2 } from "../../../modules/nf-core/bbmap/bbduk/main"

workflow BBDUK_TRIM_2PASS {

    take:
        ch_reads     // channel: [ val(meta), [ reads ] ]
        primer_file  // channel: /path/to/primer.fasta

    main:
        ch_versions      = Channel.empty()
        ch_trimmed_reads = Channel.empty()

        BBDUK_TRIM_1(ch_reads, primer_file)
        BBDUK_TRIM_2(BBDUK_TRIM_1.out.reads, primer_file)

        ch_trimmed_reads = BBDUK_TRIM_2.out.reads
        ch_versions      = ch_versions.mix(BBDUK_TRIM_2.out.versions.first())

    emit:
        trimmed_reads = ch_trimmed_reads
        versions      = ch_versions
}
