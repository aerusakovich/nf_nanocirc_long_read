process ISOCIRC {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/anrusakovich/isocirc:latest' :
        'quay.io/anrusakovich/isocirc:latest' }"

    input:
    tuple val(meta), path(fastq)
    path  fasta
    path  gtf
    path  db

    output:
    tuple val(meta), path("isocirc_output/isocirc.bed"), emit: bed
    tuple val(meta), path("isocirc_output/"),            emit: output_dir
    path  "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export LC_ALL=C
    export LANG=C
    export MPLCONFIGDIR=/tmp

    mkdir -p isocirc_output

    isocirc \\
        ${fastq} \\
        ${fasta} \\
        ${gtf} \\
        ${db} \\
        isocirc_output \\
        -t ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isocirc: \$(isocirc --version 2>&1 | head -1)
    END_VERSIONS
    """
}
