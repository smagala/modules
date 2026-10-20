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
        """
        falco ${args} \\
            -D ${prefix}_fastqc_data.txt \\
            -S ${prefix}_summary.txt \\
            -R ${prefix}_fastqc_report.html \\
            ${reads}
        """
    }
    else {
        // Falco writes to fixed output names, so a single invocation over both
        // mates overwrites R1's results with R2's. Run once per mate and derive
        // per-mate output names from the input filename instead.
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
        for read_file in ${reads}; do
            read_base=\$(basename "\$read_file")
            read_base=\${read_base%.gz}
            read_base=\${read_base%.fastq}
            read_base=\${read_base%.fq}
            falco ${args} \\
                -D "\${read_base}_fastqc_data.txt" \\
                -S "\${read_base}_summary.txt" \\
                -R "\${read_base}_fastqc_report.html" \\
                "\$read_file" &
            pids="\$pids \$!"
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
        for read_file in ${reads}; do
            read_base=\$(basename "\$read_file")
            read_base=\${read_base%.gz}
            read_base=\${read_base%.fastq}
            read_base=\${read_base%.fq}
            touch "\${read_base}_fastqc_data.txt"
            touch "\${read_base}_summary.txt"
            touch "\${read_base}_fastqc_report.html"
        done
        """
    }
}
