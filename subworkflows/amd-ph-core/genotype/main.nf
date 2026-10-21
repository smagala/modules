#!/usr/bin/env nextflow

//
// GENOTYPE: type an assembled sequence against a reference database.
//
// Generic by construction. Nothing here is pathogen-specific: the database, the thresholds and the
// meaning of a subject name are all inputs. Rabies lineage typing is this subworkflow pointed at a
// rabies typing database; flu clade calling is the same subworkflow pointed at a flu one.
//
// It also does NOT assume any particular upstream module ran. The query is an ordinary
// [ meta, fasta ] channel, so it types an ouroboros consensus, an IRMA consensus, or a FASTA a user
// dropped on disk, without knowing which.
//
// WHY THIS IS A SEPARATE STEP FROM RECRUITMENT
// --------------------------------------------
// Recruitment and typing want opposite databases. Recruitment wants FEW, well-separated references:
// reads are short, divergence is local rather than uniform, and near-identical neighbours shatter a
// sample across bins that can never be recombined, because the recruit loop only carries unmatched
// reads forward. Typing wants MANY, densely labelled references, where near-identical neighbours are
// the whole point.
//
// Measured on a real rabies library: growing the recruitment panel from 2 to 17 RABV references
// recruited 69% more reads, delivered 23% FEWER of them to the correct assembly, and manufactured a
// phantom co-infection deep enough to clear the secondary-population gate. So type AFTER assembly,
// against the dense database, where the query is a whole consensus rather than a 150 bp read and the
// local-conservation ambiguity that shatters read binning does not arise.

include { BLAST_MAKEBLASTDB } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN      } from '../../../modules/nf-core/blast/blastn/main'
include { GENOTYPE_CALL     } from '../../../modules/amd-ph-core/genotype/call/main'

workflow GENOTYPE {

    take:
        ch_query       // channel: [ val(meta), path(fasta) ]  — sequences to type
        ch_db          // channel: [ val(meta2), path(db) ]    — FASTA when build_db, else a prebuilt BLAST db dir
        build_db       // value:   boolean — build the db from ch_db, or use it as-is
        genotype_opts  // value:   map — columns, min_pident, min_qcov, tie_margin

    main:
        // A custom database is normally a FASTA someone curates, so building is the common path.
        // Accepting a prebuilt directory as well means a large shared database can be built once
        // out-of-band and reused, rather than rebuilt per run.
        if ( build_db ) {
            BLAST_MAKEBLASTDB( ch_db, [] )
            ch_blast_db = BLAST_MAKEBLASTDB.out.db
        }
        else {
            ch_blast_db = ch_db
        }

        // taxidlist / taxids / negative_tax are left empty: taxonomic filtering is a property of a
        // database that carries NCBI taxids, and a curated in-house typing set generally does not.
        // A caller that has one can filter upstream rather than have this subworkflow guess.
        BLAST_BLASTN( ch_query, ch_blast_db, [], [], false )

        GENOTYPE_CALL( BLAST_BLASTN.out.txt, genotype_opts )

    // No `versions` emit. All three components report to the `versions` TOPIC as a
    // (process, tool, version) tuple, so a pipeline collects them with one `channel.topic('versions')`
    // rather than threading a version channel out of every subworkflow it calls.
    emit:
        calls = GENOTYPE_CALL.out.calls   // channel: [ val(meta), path(*.genotype.tsv) ]
        hits  = BLAST_BLASTN.out.txt      // channel: [ val(meta), path(*.txt) ] — raw hits, for audit
}
