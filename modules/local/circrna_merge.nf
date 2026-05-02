process CIRCRNA_MERGE {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), val(tool_names), path(bed_files)
    tuple val(meta), path(pair_files)
    val   n_active

    output:
    tuple val(meta), path("${meta.id}_strict_union.bed12"),                    emit: strict_union_bed
    tuple val(meta), path("${meta.id}_strict_union_confidence.tsv"),           emit: strict_union_conf
    tuple val(meta), path("${meta.id}_strict_intersection.bed12"),             emit: strict_inter_bed
    tuple val(meta), path("${meta.id}_strict_intersection_confidence.tsv"),    emit: strict_inter_conf
    tuple val(meta), path("${meta.id}_relaxed_union.bed12"),                   emit: relaxed_union_bed
    tuple val(meta), path("${meta.id}_relaxed_union_confidence.tsv"),          emit: relaxed_union_conf
    tuple val(meta), path("${meta.id}_relaxed_intersection.bed12"),            emit: relaxed_inter_bed
    tuple val(meta), path("${meta.id}_relaxed_intersection_confidence.tsv"),   emit: relaxed_inter_conf
    path  "versions.yml",                                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def tool_args = [tool_names, bed_files.collect { it.toString() }]
        .transpose()
        .collect { name, file -> "--${name} ${file}" }
        .join(" \\\n        ")
    def pairs_arg = pair_files.collect { it.toString() }.join(" ")
    """
    # ── Step 1: BSJ-based grouping (strict + relaxed) ────────
    python3 ${projectDir}/bin/merge_circrna.py \\
        --sample  ${meta.id} \\
        ${tool_args} \\
        --outdir  .

    # ── Step 2: add isoform_confidence to all 4 TSVs ─────────
    for mode in strict relaxed; do
        for subset in union intersection; do
            python3 ${projectDir}/bin/add_isoform_confidence.py \\
                --confidence ${meta.id}_\${mode}_\${subset}_confidence.tsv \\
                --pairs      ${pairs_arg} \\
                --min_overlap ${params.circrna_isoform_overlap} \\
                --n_active   ${n_active} \\
                --output     ${meta.id}_\${mode}_\${subset}_confidence.tsv
        done
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
