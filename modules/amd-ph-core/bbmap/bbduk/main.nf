process BBMAP_BBDUK {
    tag "${meta.id}"
    label 'process_medium'

    // This is NOT a drop-in for nf-core/bbmap/bbduk, which has the same process
    // name. That one is a trimmer: one output stream, a `.log`, no `ref=`-matched
    // output and no `stats=`. This one splits the reads in two and keeps both,
    // because the matched stream is the product of a depletion step as much as
    // the depleted one is -- and it writes the stats file, which is the only
    // bbduk output MultiQC can actually parse. Installing both into one pipeline
    // would collide on the process name; install this one.
    conda "${moduleDir}/environment.yml"
    // 40.05, and not the newest build in AMDP ECR, which is
    // oamd-bio-bbtools:39.76_3ceaa9c_v3. 39.76 and 39.91 both carry a race in
    // the paired read streamer that aborts a run on two FASTQs:
    //
    //   java.lang.AssertionError: List size mismatch: 1149 vs 1148
    //       at stream.PairStreamer.nextList(PairStreamer.java:86)
    //
    // It is not cosmetic. On the pipeline's own test sample it truncated the
    // output at 2,304 of 535,442 reads and then exited non-zero; bbduk.sh runs
    // java with -ea, so the assertion is always live. 40.04 and 40.05 both
    // process the same input cleanly at one and at two threads, and produce
    // byte-identical stats files on a 38.8M-read library. ECR holds no 40.x
    // build at all and should be rebuilt at 40.05 before an AMDP config points
    // this module back at it.
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bbmap:40.05--h09cc210_0'
        : 'biocontainers/bbmap:40.05--h09cc210_0'}"

    // On MultiQC: its `bbduk` search pattern matches the string "Executing
    // jgi.BBDuk", which BBTools renamed to `bbduk.BBDukS` in 39.76, and which
    // the JVM puts on line three where the pattern allows two. It matches no
    // version this module would run, so the stderr log is not offered as a
    // MultiQC input. The `stats=` file is, under MultiQC's `bbmap/stats`
    // pattern, which content-matches the file's own header and is therefore
    // version-independent. Verified against MultiQC 1.35: the log produces "No
    // analysis results found", the stats file produces a bbmap section.
    //
    // But do not hand MultiQC the WHOLE stats file. It has one row per matched
    // reference record, and on a reference large enough to be worth using that
    // is a lot of rows: against SILVA SSU+LSU NR99, 93,601 of them, which
    // MultiQC writes verbatim into its own outputs -- 13 MB of bbmap.txt and
    // 23 MB of multiqc_data.json, for one sample. A cohort of a hundred is a
    // gigabyte of JSON. So the module writes a second, truncated copy for
    // MultiQC (`stats_mqc`) carrying the header and the top records only, and
    // leaves the complete file as the result artifact. `bbmap/stats` matches on
    // the header, so the truncated file parses identically and the summary
    // table is unchanged -- verified: 38 KB of multiqc_data against 23 MB, same
    // percent-filtered figure.

    input:
    tuple val(meta), path(reads, stageAs: "input_reads/")
    path ref

    output:
    tuple val(meta), path("*.depleted*.fastq.gz"), emit: reads
    tuple val(meta), path("*.matched*.fastq.gz") , emit: matched, optional: true
    tuple val(meta), path("*.bbduk.txt")         , emit: stats
    tuple val(meta), path("*.bbduk.summary.txt") , emit: stats_mqc
    tuple val(meta), path("*.bbduk.rpkm.txt")    , emit: rpkm
    tuple val(meta), path("*.bbduk.refstats.txt"), emit: refstats
    tuple val("${task.process}"), val('bbduk'), eval("bbduk.sh --version 2>&1 | sed -n 's/^BBTools version //p' | head -n 1"), topic: versions, emit: versions_bbduk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Sorted so mate 1 is always first: the glob that stages the reads makes no
    // ordering promise, and passing them the wrong way round would swap the
    // output mates silently. Matches how DEACON_FILTER and HOSTESS do it.
    def sorted_reads = meta.single_end ? [reads].flatten() : reads.sort { read -> read.simpleName }
    def in_cmd  = meta.single_end
        ? "in=${sorted_reads[0]}"
        : "in=${sorted_reads[0]} in2=${sorted_reads[1]}"
    // `out` is what did NOT match the reference and `outm` is what did. For this
    // pipeline the depleted stream is the product and the matched stream is the
    // input to rRNA characterisation, so both are kept.
    def out_cmd = meta.single_end
        ? "out=${prefix}.depleted.fastq.gz outm=${prefix}.matched.fastq.gz"
        : "out=${prefix}.depleted_1.fastq.gz out2=${prefix}.depleted_2.fastq.gz " +
          "outm=${prefix}.matched_1.fastq.gz outm2=${prefix}.matched_2.fastq.gz"
    // BBTools sizes its heap from the machine, not from the task, and will
    // happily ask for more than the executor granted and be killed for it. 85%
    // rather than all of it: the JVM's non-heap footprint (metaspace, thread
    // stacks, direct buffers) sits outside -Xmx, so handing it the whole cgroup
    // limit is zero-slack by construction. Measured peak RSS for bbduk is
    // -Xmx plus ~0.1 GB, so 85% is comfortable rather than tight. Floored at
    // 1g so a small test task still gets a legal value.
    def heap_gb = task.memory ? Math.max(1L, (task.memory.toGiga() * 0.85) as long) : 0L
    def mem = task.memory ? "-Xmx${heap_gb}g" : ''
    // Rows kept in the MultiQC-facing copy of the stats file. 20 is MultiQC's
    // own default for how many bars it draws in that plot, so more would be
    // carried and not shown. Not exposed as an ext key: nf-core lint allows
    // only args/args2/args3/prefix/when/suffix, and nothing needs to tune it.
    def mqc_rows = 20
    """
    bbduk.sh \\
        ${mem} \\
        ${in_cmd} \\
        ${out_cmd} \\
        ref=${ref} \\
        stats=${prefix}.bbduk.txt \\
        rpkm=${prefix}.bbduk.rpkm.txt \\
        refstats=${prefix}.bbduk.refstats.txt \\
        threads=${task.cpus} \\
        ${args}

    # Header (4 lines) plus the top `mqc_rows` records. The records are already
    # in descending read order, so this is the head of the distribution, not a
    # sample of it.
    awk 'NR <= 4 || n++ < ${mqc_rows}' ${prefix}.bbduk.txt > ${prefix}.bbduk.summary.txt
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def outputs = meta.single_end
        ? ["${prefix}.depleted.fastq.gz", "${prefix}.matched.fastq.gz"]
        : ["${prefix}.depleted_1.fastq.gz", "${prefix}.depleted_2.fastq.gz",
           "${prefix}.matched_1.fastq.gz", "${prefix}.matched_2.fastq.gz"]
    """
    for f in ${outputs.join(' ')}; do
        echo "" | gzip > \$f
    done

    printf '%s\\n' \\
        '#File\tstub' \\
        '#Total\t0' \\
        '#Matched\t0\t0.00000%' \\
        '#Name\tReads\tReadsPct' > ${prefix}.bbduk.txt

    cp ${prefix}.bbduk.txt ${prefix}.bbduk.summary.txt

    printf '%s\\n' \\
        '#File\tstub' \\
        '#Reads\t0' \\
        '#Mapped\t0' \\
        '#RefSequences\t0' \\
        '#Name\tLength\tBases\tCoverage\tReads\tRPKM' > ${prefix}.bbduk.rpkm.txt

    printf '%s\\n' \\
        '#File\tstub' \\
        '#Reads\t0' \\
        '#Mapped\t0' \\
        '#References\t0' \\
        '#Name\tLength\tScaffolds\tBases\tCoverage\tReads\tRPKM' > ${prefix}.bbduk.refstats.txt
    """
}
