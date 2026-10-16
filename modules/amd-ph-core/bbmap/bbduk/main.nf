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
    // 40.04, and not the newest build in AMDP ECR, which is
    // oamd-bio-bbtools:39.76_3ceaa9c_v3. 39.76 and 39.91 both carry a race in
    // the paired read streamer that aborts a run on two FASTQs:
    //
    //   java.lang.AssertionError: List size mismatch: 1149 vs 1148
    //       at stream.PairStreamer.nextList(PairStreamer.java:86)
    //
    // It is not cosmetic. On the pipeline's own test sample it truncated the
    // output at 2,304 of 535,442 reads and then exited non-zero; bbduk.sh runs
    // java with -ea, so the assertion is always live. 40.04 processes the same
    // input cleanly at one and at two threads. ECR should be rebuilt at 40.04
    // before an AMDP config points this module back at it.
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bbmap:40.04--h09cc210_0'
        : 'biocontainers/bbmap:40.04--h09cc210_0'}"

    // On MultiQC: its `bbduk` search pattern matches the string "Executing
    // jgi.BBDuk", which BBTools renamed to `bbduk.BBDukS` in 39.76, and which
    // the JVM puts on line three where the pattern allows two. It matches no
    // version this module would run, so the stderr log is not offered as a
    // MultiQC input. The `stats=` file is, under MultiQC's `bbmap/stats`
    // pattern, which content-matches the file's own header and is therefore
    // version-independent. Verified against MultiQC 1.35: the log produces "No
    // analysis results found", the stats file produces a bbmap section.

    input:
    tuple val(meta), path(reads, stageAs: "input_reads/")
    path ref

    output:
    tuple val(meta), path("*.depleted*.fastq.gz"), emit: reads
    tuple val(meta), path("*.matched*.fastq.gz") , emit: matched, optional: true
    tuple val(meta), path("*.bbduk.txt")         , emit: stats
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
    // happily ask for more than the executor granted and be killed for it.
    def mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : ''
    """
    bbduk.sh \\
        ${mem} \\
        ${in_cmd} \\
        ${out_cmd} \\
        ref=${ref} \\
        stats=${prefix}.bbduk.txt \\
        threads=${task.cpus} \\
        ${args}
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
    """
}
