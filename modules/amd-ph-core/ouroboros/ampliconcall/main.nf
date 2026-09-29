// OUROBOROS_AMPLICONCALL — reference-anchored, coverage-masked output for amplicon and patchy-coverage
// samples.
//
// The iterative assembler force-extends a contiguous consensus across thin (1-4x) amplicon coverage and
// splices out true-zero gaps, which loses genome coordinates. For amplicon and RACE libraries — uneven
// coverage, dropout regions — that is misleading. This instead maps the recruited reads to the gene's
// FULL panel reference and emits, anchored to genome coordinates:
//   <gene>.amplicon.bam(.bai)      reads vs the full reference: the island view, at true depth
//   <gene>.amplicon.consensus.fa   full-length consensus: a base where depth >= min_depth, else N
//   <gene>.amplicon.islands.bed    intervals with depth >= min_depth (the covered islands)
//
// Honest by construction: it never bridges a gap it has no reads for, and never collapses one.
//
// The reads are INFLATED here rather than used as patterns, because a pattern count understates depth
// and every threshold below is a per-read depth.

process OUROBOROS_AMPLICONCALL {
    tag "${meta.id}_${gene}"
    label 'process_medium'

    container 'oamd-bio-python:3.12-orb'

    input:
    tuple val(meta), val(gene), path(reads), path(deflation_index), path(refs)
    val min_depth

    output:
    tuple val(meta), val(gene), path("${gene}.amplicon.bam"), path("${gene}.amplicon.bam.bai"), emit: bam
    tuple val(meta), val(gene), path("${gene}.amplicon.consensus.fa")                         , emit: consensus
    tuple val(meta), val(gene), path("${gene}.amplicon.islands.bed")                          , emit: bed
    tuple val(meta), val(gene), path("${gene}.ref")                                           , emit: ref
    path "versions.yml"                                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Inflate the recruited read patterns to real reads, so depth below is TRUE per-read depth.
    cat ${reads} > ${gene}_all.fa
    irma-core xflate --inflate ${deflation_index} ${gene}_all.fa > ${gene}.fastq

    # The FULL-length panel record for this gene: the genome-coordinate anchor. The record ID is the
    # gene/bin label, which is what makes reference-framed reporting possible downstream (ADR-0014).
    # Match on the BIN, read from the defline's `bin=` tag. Record ids are accessions, so comparing
    # an id against a gene name never matches; a record with no `bin=` tag is its own bin, which is
    # what the refined per-gene references the gather loop writes need. This previously compared a
    # `{Sxx}`-suffixed id against the gene, emitted an EMPTY .ref and still exited 0 — the failure
    # only surfaced downstream in STITCH.
    awk -v g="${gene}" '/^>/{ n=split(substr(\$0,2),a," "); b=a[1];
                              for(i=2;i<=n;i++) if(a[i] ~ /^bin=/) b=substr(a[i],5);
                              p=(b==g) } p' ${refs} > ${gene}.ref
    if [ ! -s "${gene}.ref" ]; then
        echo "ERROR: no panel record for bin ${gene} in ${refs}" >&2
        exit 1
    fi

    # Map the inflated reads to the full reference, sort and index: the island BAM.
    rammap -a -x sr -t ${task.cpus} ${gene}.ref ${gene}.fastq 2>/dev/null \\
        | samtools sort -@ ${task.cpus} -o ${gene}.amplicon.bam -
    samtools index ${gene}.amplicon.bam

    # Reference-anchored, gap-masked consensus: a real base where depth >= min, N elsewhere.
    samtools consensus -a -f fasta -d ${min_depth} -o ${gene}.amplicon.consensus.fa ${gene}.amplicon.bam

    # Covered islands (depth >= min) as a BED.
    samtools depth -a ${gene}.amplicon.bam | awk -v d=${min_depth} 'BEGIN{OFS="\\t"}
        { cov=(\$3>=d)?1:0
          if(cov && !inq){st=\$2-1; inq=1}
          if(!cov && inq){print \$1, st, \$2-1; inq=0} }
        END{ if(inq) print \$1, st, \$2 }' > ${gene}.amplicon.islands.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma_core: \$(irma-core --version 2>&1 | grep -oE '[0-9][^ ]*' | head -1)
        rammap: \$(rammap --version 2>/dev/null | grep -oE '[0-9][^ ]*' | head -1)
        samtools: \$(samtools --version 2>/dev/null | head -1 | grep -oE '[0-9.]+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${gene}.amplicon.bam ${gene}.amplicon.bam.bai ${gene}.amplicon.consensus.fa \\
          ${gene}.amplicon.islands.bed ${gene}.ref

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma_core: \$(irma-core --version 2>&1 | grep -oE '[0-9][^ ]*' | head -1)
        rammap: \$(rammap --version 2>/dev/null | grep -oE '[0-9][^ ]*' | head -1)
        samtools: \$(samtools --version 2>/dev/null | head -1 | grep -oE '[0-9.]+' | head -1)
    END_VERSIONS
    """
}
