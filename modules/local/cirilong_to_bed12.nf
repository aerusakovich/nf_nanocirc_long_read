process CIRILONG_TO_BED12 {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(info)      // CIRI-long.info file

    output:
    tuple val(meta), path("${meta.id}_cirilong.bed12"), emit: bed12
    path  "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${projectDir}/bin/cirilong_to_bed12.py \\
        ${info} \\
        ${meta.id}_cirilong.bed12

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
