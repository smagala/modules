// OUROBOROS_GATHERREPORTS — collect the recruit loop's terminal reporting files for one sample:
//   chimeric.fa            round-1 reads dropped as chimeric
//   nomatch.fa             reads still unrecruited after the final round
//   sorted_read_stats.txt  per-gene read-pattern and read counts, every round, tab-separated
//
// Inspection outputs only; nothing downstream consumes them. nomatch.fa in particular is what lets a
// reader distinguish "not recruited" from "not sequenced", which matters because every depth this
// pipeline reports is recruited depth.
//
// No versions.yml: this process only copies and concatenates files.

process OUROBOROS_GATHERREPORTS {
    tag "${meta.id}"
    label 'process_low'

    container 'oamd-bio-ouroboros:1.0.0'

    input:
    tuple val(meta), path(chimeric), path(nomatch), path(sort_stats, stageAs: 'sort_stats_*.txt')

    output:
    tuple val(meta), path("chimeric.fa"), path("nomatch.fa"), path("sorted_read_stats.txt"), emit: reports

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cp ${chimeric} chimeric.fa
    cp ${nomatch} nomatch.fa
    printf 'Gene\\tRead Patterns\\tRead Count\\n' > sorted_read_stats.txt
    cat ${sort_stats} >> sorted_read_stats.txt
    """

    stub:
    """
    touch chimeric.fa nomatch.fa
    printf 'Gene\\tRead Patterns\\tRead Count\\n' > sorted_read_stats.txt
    """
}
