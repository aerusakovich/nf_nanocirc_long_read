process CIRCNICK_TO_BED12 {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(annotated)    // circRNA_candidates.annotated.txt
                                        // (may be lifted if liftover was run)
    tuple val(meta), path(exon_usage)   // circ_circRNA_exon_usage_length_of_exons.txt
                                        // (may be lifted if liftover was run)
    tuple val(meta), path(intron_cov)   // introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.intronCov.bed
                                        // used to recover exon structure for circRNAs
                                        // missing from or with empty coords in exon_usage

    output:
    tuple val(meta), path("${meta.id}_circnick.bed12"), emit: bed12
    path  "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${projectDir}/bin/circnick_to_bed12.py \\
        ${annotated} \\
        ${exon_usage} \\
        ${intron_cov} \\
        ${meta.id}_circnick.bed12

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
