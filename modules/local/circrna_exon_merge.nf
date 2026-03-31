process CIRCRNA_EXON_MERGE {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(pair_files)
    val   n_active

    output:
    tuple val(meta), path("${meta.id}_exon_union.bed12"),                emit: exon_union_bed
    tuple val(meta), path("${meta.id}_exon_union_confidence.tsv"),       emit: exon_union_conf
    tuple val(meta), path("${meta.id}_exon_intersection.bed12"),         emit: exon_inter_bed
    tuple val(meta), path("${meta.id}_exon_intersection_confidence.tsv"),emit: exon_inter_conf
    path  "versions.yml",                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def pairs_arg = pair_files.collect { it.toString() }.join(" \\\n        ")
    """
    python3 ${projectDir}/bin/merge_exon_based.py \\
        --sample      ${meta.id} \\
        --pairs       ${pairs_arg} \\
        --min_overlap ${params.circrna_isoform_overlap} \\
        --n_active    ${n_active} \\
        --outdir      .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
