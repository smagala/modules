// OUROBOROS_GENOMEREPORT — per-sample genomes / contamination screen.
//
// For every assembled genome of a sample, aligns its consensus to ITS OWN panel record and tabulates
// read count, consensus length, %called, identity to that record, and aligned fraction. That last pair
// is the realness annotation: it separates a genuine co-infection or contaminant from a conserved-region
// cross-map echo of the primary, which otherwise looks like a second genome.
//
// This is a CROSS-BIN operation — it compares every genome of a sample against every other — so it
// cannot live inside the per-bin assemble scatter. It runs once per sample, after all bins are in.

process OUROBOROS_GENOMEREPORT {
    tag "${meta.id}"
    label 'process_low'

    container 'oamd-bio-python:3.12-orb'

    input:
    tuple val(meta), path(consensuses), path(panel_ref), path(sort_stats)

    output:
    tuple val(meta), path("*.genomes.tsv"), emit: report
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ss = (sort_stats && sort_stats.name != 'NO_FILE') ? "--sort-stats ${sort_stats}" : ''
    """
    # The ref-pinned consensuses are the inputs — everything except the panel reference itself. Take the
    # target name from each FASTA header rather than the filename, so this is robust to both the
    # <gene>.stitched.fa and <gene>.fa naming, then align each to its OWN panel record.
    for c in ${consensuses}; do
        [ "\$c" = "${panel_ref}" ] && continue
        tg=\$(head -1 "\$c" | sed 's/^>//; s/[ \\t].*//')
        # Match on the record id with its `{Sxx}` variant suffix STRIPPED. The consensus header is a
        # bin/gene name (`..._khatangaense_L`) while the panel writes competing variants as
        # `..._khatangaense_L{S01}`, so an exact-or-prefix compare matched nothing, wrote an empty
        # ref.fa, and left rammap emitting a headerless SAM. That is why identity_to_ref and
        # aligned_frac have been NA/0.00 on every row this module has ever produced — and why
        # pct_called silently fell back to the unbounded called/ref_len denominator.
        #
        # EVERY record of the bin, not just the first. A bin holds competing variants precisely
        # because one of them fits the sample and the others do not, and the assembly was built
        # against whichever won — so measuring identity against the bin's first record measures the
        # wrong sequence. It reported khatangaense_L at 81.7% (vs NC_055196) when the assembly is
        # 100.0% to MK352484, and West Nile at 78.8% when it is 99.8% to EF429197. That made 25 of
        # 85 rows look like sub-90% "nearest neighbour" calls when only 5 are. Emitting all records
        # lets rammap pick the best and identity_to_ref resolve per alignment.
        awk -v n="\$tg" '/^>/{ k=split(substr(\$0,2),a," "); b=a[1];
                               for(i=2;i<=k;i++) if(a[i] ~ /^bin=/) b=substr(a[i],5);
                               p=(b==n) } p' ${panel_ref} > "\$tg.ref.fa"
        rammap -a --eqx "\$tg.ref.fa" "\$c" 2>/dev/null > "\$tg.sam" || : > "\$tg.sam"
    done

    CONS=\$(for c in ${consensuses}; do [ "\$c" = "${panel_ref}" ] || echo "\$c"; done)
    genome_report.py \$CONS \\
        --sample ${meta.id} \\
        --ref-fasta ${panel_ref} \\
        --sam-dir . \\
        ${ss} > ${meta.id}.genomes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rammap: \$(rammap --version 2>/dev/null | grep -oE '[0-9][^ ]*' | head -1 || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.genomes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rammap: \$(rammap --version 2>/dev/null | grep -oE '[0-9][^ ]*' | head -1 || echo unknown)
    END_VERSIONS
    """
}
