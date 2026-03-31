process FASTQ_RENAME {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: fastq

    script:
    """
    ln -s ${fastq} ${meta.id}.fastq.gz
    """
}
