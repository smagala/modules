process FALCO {
    tag "${meta.id}"
    // process_low rather than process_single: paired-end input runs one falco
    // per mate concurrently, so the task wants two cores. Single-end input
    // leaves the second idle, which is the cheaper mistake — the alternative is
    // a one-core reservation that paired-end input immediately doubles.
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/falco:1.2.5--h077b44d_0'
        : 'biocontainers/falco:1.2.5--h077b44d_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc_report.html"), emit: html
    tuple val(meta), path("*_fastqc_data.txt")   , emit: txt
    tuple val(meta), path("*_summary.txt")       , emit: summary
    tuple val("${task.process}"), val('falco'), eval("falco --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_falco

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // --threads is deliberately not passed. falco accepts -t/--threads only as
    // a FastQC compatibility stub — its own help marks it "[NOT YET IMPLEMENTED
    // IN FALCO]" — so it has never done anything here. It is dropped rather
    // than left in place because the meaning it would have if implemented is
    // FastQC's: the number of *files* processed at once, which is what the loop
    // below already does. Handing task.cpus to each of two concurrent
    // invocations would then ask for twice the allocation.
    if (reads.toList().size() == 1) {
        // Linked to the prefix for the same reason as the paired branch below:
        // falco records the name of the file it was handed as `Filename`, and
        // MultiQC's FastQC parser takes its sample name from that rather than
        // from the report's path.
        """
        case "${reads}" in
            *.gz) ext=".fastq.gz" ;;
            *)    ext=".fastq" ;;
        esac
        ln -s "\$(readlink -f "${reads}")" "${prefix}\${ext}"

        falco ${args} \\
            -D ${prefix}_fastqc_data.txt \\
            -S ${prefix}_summary.txt \\
            -R ${prefix}_fastqc_report.html \\
            "${prefix}\${ext}"
        """
    }
    else {
        // Falco writes to fixed output names, so a single invocation over both
        // mates overwrites R1's results with R2's. Run once per mate, naming
        // everything `<prefix>_R1` and `<prefix>_R2` from the mate's position
        // in the staged list rather than from its filename, so the same sample
        // reports under the same name whatever the sequencer or an upstream
        // step called the file.
        //
        // Falco is given a LINK named for the mate rather than the staged file
        // itself, and this is the part that actually matters. Falco records the
        // name of the file it was handed into the report, as `Filename`, and
        // MultiQC's FastQC parser takes its sample name from that embedded
        // value -- not from the report's path. Renaming only the outputs would
        // leave the reports correctly named on disk and still grouped under the
        // input's name in MultiQC, which is the bug this is meant to fix. The
        // extension is carried over so falco still detects compression.
        //
        // The mates run concurrently rather than one after the other. falco is
        // single-threaded, so this is the only parallelism available, and it
        // roughly halves wall time on paired-end input. A second concurrent
        // invocation costs little memory: falco's footprint is set by its fixed
        // per-file tables rather than by the size of the input.
        //
        // Each PID is waited on individually so its status is checked. A bare
        // `wait` returns zero however the children exited, which would leave a
        // missing or truncated report behind a task that reported success.
        """
        pids=''
        mate=1
        for read_file in ${reads}; do
            read_base="${prefix}_R\${mate}"
            case "\$read_file" in
                *.gz) ext=".fastq.gz" ;;
                *)    ext=".fastq" ;;
            esac
            ln -s "\$(readlink -f "\$read_file")" "\${read_base}\${ext}"
            falco ${args} \\
                -D "\${read_base}_fastqc_data.txt" \\
                -S "\${read_base}_summary.txt" \\
                -R "\${read_base}_fastqc_report.html" \\
                "\${read_base}\${ext}" &
            pids="\$pids \$!"
            mate=\$((mate + 1))
        done

        for pid in \$pids; do
            wait "\$pid"
        done
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (reads.toList().size() == 1) {
        """
        touch ${prefix}_fastqc_data.txt
        touch ${prefix}_summary.txt
        touch ${prefix}_fastqc_report.html
        """
    }
    else {
        """
        for mate in 1 2; do
            read_base="${prefix}_R\${mate}"
            touch "\${read_base}_fastqc_data.txt"
            touch "\${read_base}_summary.txt"
            touch "\${read_base}_fastqc_report.html"
        done
        """
    }
}
