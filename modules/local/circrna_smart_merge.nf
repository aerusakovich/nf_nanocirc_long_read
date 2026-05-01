process CIRCRNA_SMART_MERGE {
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
    tuple val(meta), path("${meta.id}_discovery.bed12"),                          emit: hybrid_bed
    tuple val(meta), path("${meta.id}_discovery_confidence.tsv"),                 emit: hybrid_conf
    tuple val(meta), path("${meta.id}_smart_consensus.bed12"),                    emit: consensus_bed
    tuple val(meta), path("${meta.id}_smart_consensus_confidence.tsv"),           emit: consensus_conf
    tuple val(meta), path("${meta.id}_smart_consensus_xstruct.bed12"),            emit: consensus_xstruct_bed
    tuple val(meta), path("${meta.id}_smart_consensus_xstruct_confidence.tsv"),   emit: consensus_xstruct_conf
    tuple val(meta), path("${meta.id}_smart_priority.bed12"),                     emit: priority_bed
    tuple val(meta), path("${meta.id}_smart_priority_confidence.tsv"),            emit: priority_conf
    path  "versions.yml",                                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def names_arg = tool_names.join(' ')
    def files_arg = bed_files.collect { it.toString() }.join(' ')
    def pairs_arg = pair_files.collect { it.toString() }.join(' ')
    """
    python3 ${projectDir}/bin/smart_merge.py \\
        --sample           ${meta.id} \\
        --tool_names       ${names_arg} \\
        --bed_files        ${files_arg} \\
        --tolerance        ${params.circrna_bsj_tolerance} \\
        --struct_tolerance ${params.circrna_bsj_tolerance} \\
        --n_active         ${n_active} \\
        --outdir           .

    # Add Low/Medium/High confidence scoring to all four mode TSVs
    for mode in consensus consensus_xstruct consensus_hybrid priority; do
        tsv="${meta.id}_smart_\${mode}_confidence.tsv"
        python3 ${projectDir}/bin/add_isoform_confidence.py \\
            --confidence  "\${tsv}" \\
            --pairs       ${pairs_arg} \\
            --min_overlap ${params.circrna_isoform_overlap} \\
            --n_active    ${n_active} \\
            --strip_isoform_suffix \\
            --output      "\${tsv}"
    done

    # Rename hybrid outputs to user-facing discovery mode names
    mv ${meta.id}_smart_consensus_hybrid.bed12           ${meta.id}_discovery.bed12
    mv ${meta.id}_smart_consensus_hybrid_confidence.tsv  ${meta.id}_discovery_confidence.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
