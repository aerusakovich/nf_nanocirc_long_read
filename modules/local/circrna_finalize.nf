process CIRCRNA_FINALIZE {
    tag "${meta.id}_${meta.category}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pybedtools:0.12.0--py39h475c85d_0' :
        'quay.io/biocontainers/pybedtools:0.12.0--py39h475c85d_0' }"

    input:
    tuple val(meta), path(annotated_tsv), path(iso_expr), path(ciri_expr), path(fl_expr), path(nick_expr)
    path gene_bed
    path exon_bed

    output:
    tuple val(meta), path("*_clean.tsv"), emit: clean
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}_${meta.category}"
    def iso_arg  = iso_expr.name  != 'NO_FILE' ? "--iso_expr    ${iso_expr}"   : ''
    def ciri_arg = ciri_expr.name != 'NO_FILE' ? "--ciri_expr   ${ciri_expr}"  : ''
    def fl_arg   = fl_expr.name   != 'NO_FILE' ? "--circfl_expr ${fl_expr}"    : ''
    def nick_arg = nick_expr.name != 'NO_FILE' ? "--nick_expr   ${nick_expr}"  : ''
    """
    circrna_clean.py \\
        --annotated_tsv ${annotated_tsv} \\
        --gene_bed      ${gene_bed} \\
        --exon_bed      ${exon_bed} \\
        --prefix        ${prefix} \\
        --bsj_tol       ${params.circrna_bsj_tolerance} \\
        ${iso_arg} \\
        ${ciri_arg} \\
        ${fl_arg} \\
        ${nick_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        bedtools: \$(bedtools --version | head -1 | sed 's/bedtools v//')
    END_VERSIONS
    """
}
