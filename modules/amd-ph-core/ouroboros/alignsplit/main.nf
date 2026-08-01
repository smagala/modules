// OUROBOROS_ALIGNSPLIT — partition the per-gene (already-inflated) reads into chunk FASTQs.
//
// First stage of the alignment scatter-gather (ALIGNSPLIT -> ALIGNCHUNK -> ALIGNMERGE). ALIGNCHUNK
// threads over a chunk's reads (multicore), so a sample at or below `chunk_reads` stays ONE fat
// multicore task; only beyond it do we split into even ~`chunk_reads` chunks that fan out across cloud
// instances, where the compute saved outweighs per-instance startup. RAMMAP threads internally and is
// always one chunk. Emits the chunk count so the ALIGNMERGE gather can be bounded (groupKey) under the
// assemble recursion — an unbounded groupTuple hangs a Nextflow recursion.

process OUROBOROS_ALIGNSPLIT {
    tag "${meta.id}_${gene}_r${round}"
    label 'process_low'

    container 'oamd-bio-ouroboros:1.0.0'

    input:
    tuple val(meta), val(gene), val(round), path(reads), path(ref)
    val aligner
    val chunk_reads

    output:
    tuple val(meta), val(gene), val(round), path(ref), path("chunk.*.fastq"), emit: chunks
    tuple val(meta), val(gene), val(round), env('NCHUNK')                   , emit: nchunk
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat ${reads} > all.fastq
    nrec=\$(( \$(wc -l < all.fastq) / 4 ))

    per=${chunk_reads}
    [ "\$per" -lt 1 ] && per=500000
    if [ "${aligner}" = "RAMMAP" ] || [ "\$nrec" -le "\$per" ]; then
        # BYPASS the scatter at or below chunk_reads — one fat MULTICORE task (ALIGNCHUNK threads over
        # the chunk), no fan-out, no per-instance startup. Also always for RAMMAP (threads internally).
        # Beyond chunk_reads we scatter, because a single task's compute then exceeds cloud startup
        # overhead by enough that spreading even chunks across instances wins wall time.
        mv all.fastq chunk.0000.fastq
        NCHUNK=1
    else
        # EVEN split: ceil(nrec/per) chunks of ~equal size, not full chunks plus a small remainder.
        # Reads are per-read independent, so any partition merges identically.
        nch=\$(( (nrec + per - 1) / per ))
        per_even=\$(( (nrec + nch - 1) / nch ))
        split -d -a 4 -l \$(( per_even * 4 )) all.fastq chunk. --additional-suffix=.fastq
        rm -f all.fastq
        NCHUNK=\$(ls chunk.*.fastq | wc -l)
        [ "\$NCHUNK" -lt 1 ] && { : > chunk.0000.fastq; NCHUNK=1; }
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(split --version | head -1 | grep -oE '[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch chunk.0000.fastq
    NCHUNK=1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(split --version | head -1 | grep -oE '[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """
}
