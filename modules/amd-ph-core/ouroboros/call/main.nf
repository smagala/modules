// OUROBOROS_CALL — terminal variant calling and consensus from the converged reference and its SAM.
//
// The non-iterated back half of assembly: BAM; on the paired path irma-core merge-sam, then
// var_call_stats and call; phasing when there are more than two variants; amended consensus; VCF.
// One process per gene, so segments call in parallel.
//
// NOTE (ADR-0014): the VCF produced here is in the ASSEMBLED CONSENSUS frame, so it reports only
// within-sample (minority / subclonal) variation. Fixed differences from the panel record are absent
// by construction, because reference refinement moved the reference onto them before calling.
// Reference-framed reporting is a separate, additive step.

process OUROBOROS_CALL {
    tag "${meta.id}_${gene}"
    label 'process_medium'

    container 'oamd-bio-python:3.12-orb'

    input:
    tuple val(meta), val(gene), path(ref), path(sam)
    val min_variant_count
    val min_variant_freq
    val min_ins_freq
    val min_del_freq
    val min_variant_qual
    val min_column_coverage
    val min_confidence
    val sig_level
    val auto_freq
    val min_ambig_freq
    val seg_numbers
    val sort_groups
    val caller
    val refchg
    val sor
    val amend_indels
    val refchg_min_strand_frac

    output:
    tuple val(meta), val(gene), path("${gene}.bam")              , emit: bam
    tuple val(meta), val(gene), path("${gene}.bam.bai")          , emit: bai
    tuple val(meta), val(gene), path("${gene}.fasta")            , emit: consensus
    tuple val(meta), val(gene), path("${gene}-variants.txt")     , emit: variants
    tuple val(meta), val(gene), path("${gene}-allAlleles.txt")   , emit: alleles
    tuple val(meta), val(gene), path("${gene}-coverage.txt")     , emit: coverage
    tuple val(meta), val(gene), path("amended_consensus/*.fa")   , emit: amended    , optional: true
    tuple val(meta), val(gene), env('AMENDED_LEN')               , emit: amended_len
    tuple val(meta), val(gene), path("${gene}-insertions.txt")   , emit: insertions , optional: true
    tuple val(meta), val(gene), path("${gene}-deletions.txt")    , emit: deletions  , optional: true
    tuple val(meta), val(gene), path("${gene}.vcf")              , emit: vcf        , optional: true
    tuple val(meta), val(gene), path("${gene}-pairingStats.txt") , emit: pairing    , optional: true
    tuple val(meta), val(gene), path("${gene}*.sqm")             , emit: sqm        , optional: true
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def paired    = meta.single_end ? 0 : 1
    def seg_opt   = seg_numbers ? "--seg ${seg_numbers}" : ''
    // Header numbering applies only when the panel is NOT sort-grouped, i.e. non-segmented.
    def h_opt     = sort_groups ? '' : '--fa-header-suffix'
    // Opt-in indel folding: fold called indels into the amended consensus so it reflects real length
    // changes (N still means missing). Off by default, because IRMA routes indels to the VCF and keeps
    // the amended consensus length-locked; folding here would diverge the NM=0 parity baseline.
    // Optional consensus depth floor. Applies to call.py ONLY -- CALL_OPTS is shared with
    // generate_vcf.py, which does not know this flag. Defaults off: existing callers unchanged.
    def mcd = ((task.ext.min_consensus_depth ?: 0) as int) > 0
                  ? "--min-consensus-depth ${task.ext.min_consensus_depth} " : ''
    def indel_opt = amend_indels ? "--deletion-file ${gene}-deletions.txt --insertion-file ${gene}-insertions.txt" : ''
    """
    # CALL_OPTS holds the thresholds shared by call.py AND generate_vcf.py. The caller-specific flags
    # stay out of it: generate_vcf.py's argparse rejects them.
    CALL_OPTS="--min-count ${min_variant_count} --min-freq ${min_variant_freq} --min-insertion-freq ${min_ins_freq} --min-deletion-freq ${min_del_freq} --min-quality ${min_variant_qual} --min-total-col-coverage ${min_column_coverage} --conf-not-mac-err ${min_confidence} --sig-level ${sig_level}"
    ${auto_freq == 1 ? 'CALL_OPTS="\$CALL_OPTS --auto-min-freq"' : 'true'}
    CALLER_OPTS="--caller ${caller} --sor ${sor} --refchg-min-strand-frac ${refchg_min_strand_frac}"

    # === SAM -> BAM ===
    samtools view -bS ${sam} > ${gene}.bam 2>/dev/null
    samtools sort ${gene}.bam -o ${gene}.bam 2>/dev/null
    samtools index ${gene}.bam 2>/dev/null

    # === VARIANT CALLING ===
    if [ "${paired}" -eq 1 ]; then
        # Merge mate pairs and pileup on the merged SAM; the BAM above stays on the raw SAM.
        irma-core merge-sam --store-stats ${ref} ${sam} M-${gene}
        var_call_stats.py ${ref} M-${gene}.sam V-${gene}
        get_pairing_stats.py M-${gene}.stats > ${gene}-pairingStats.txt
        call.py --print-all-sites --no-gap-allele ${mcd}\$CALL_OPTS \$CALLER_OPTS --paired-error ${gene}-pairingStats.txt ${ref} ${gene} V-${gene}.pileup.msgpack
    else
        var_call_stats.py ${ref} ${sam} V-${gene}
        call.py --print-all-sites --no-gap-allele ${mcd}\$CALL_OPTS \$CALLER_OPTS ${ref} ${gene} V-${gene}.pileup.msgpack
    fi

    # === PHASING === (only worthwhile with more than two variants)
    n=\$(wc -l < ${gene}-variants.txt || echo 0)
    if [ "\$n" -gt 2 ]; then
        phase.py ${gene} ${gene}-vars.msgpack ${gene}-pats.msgpack --array-size 1 --index 1
        for metric in EXPENRD JACCARD MUTUALD NJOINTP; do
            [ -f ${gene}-\${metric}.sqm ] && complete_matrix.py ${gene}-\${metric}.sqm
        done
    fi

    # === CONSENSUS ===
    mkdir -p amended_consensus
    amend_consensus.py --prefix amended_consensus --count ${min_variant_count} --freq ${min_ambig_freq} --refchg ${refchg} ${indel_opt} ${seg_opt} ${h_opt} ${gene}.fasta ${gene}-variants.txt

    # === VCF ===
    if [ -f ${gene}-insertions.txt ] && [ -f ${gene}-deletions.txt ]; then
        generate_vcf.py \$CALL_OPTS ${gene}.fasta ${gene}-allAlleles.txt ${gene}-insertions.txt ${gene}-deletions.txt --out ${gene}.vcf
    fi

    # Amended sequence length (0 if none). The stitch gate downstream relabels amended -> stitched only
    # when the amended consensus already spans the full panel reference, so it needs this as a value.
    if ls amended_consensus/*.fa >/dev/null 2>&1; then
        AMENDED_LEN=\$(grep -hv '^>' amended_consensus/*.fa | tr -d '\\n\\r' | wc -c)
    else
        AMENDED_LEN=0
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma_core: \$(irma-core --version 2>&1 | grep -oE '[0-9][^ ]*' | head -1)
        samtools: \$(samtools --version 2>/dev/null | head -1 | grep -oE '[0-9.]+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p amended_consensus
    touch ${gene}.bam ${gene}.bam.bai ${gene}.fasta ${gene}-variants.txt ${gene}-allAlleles.txt ${gene}-coverage.txt
    AMENDED_LEN=0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma_core: \$(irma-core --version 2>&1 | grep -oE '[0-9][^ ]*' | head -1)
        samtools: \$(samtools --version 2>/dev/null | head -1 | grep -oE '[0-9.]+' | head -1)
    END_VERSIONS
    """
}
