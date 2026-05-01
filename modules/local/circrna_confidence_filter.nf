process CIRCRNA_CONFIDENCE_FILTER {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(bed12), path(tsv)

    output:
    tuple val(meta), path("*.bed12"), emit: bed
    tuple val(meta), path("*_confidence.tsv"), emit: conf
    path  "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}_${meta.category}"
    def mode    = task.ext.args   ?: 'no_low'
    """
    python3 ${projectDir}/bin/filter_confidence.py \\
        --bed    ${bed12} \\
        --tsv    ${tsv} \\
        --mode   ${mode} \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
