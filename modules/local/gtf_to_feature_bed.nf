process GTF_TO_FEATURE_BED {
    tag "gtf"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    path gtf

    output:
    path "gene.bed",     emit: gene_bed
    path "exon.bed",     emit: exon_bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk 'BEGIN{OFS="\\t"} !/^#/ && \$3=="gene" {
        n="."
        for(i=9;i<=NF;i++){if(\$i=="gene_id"){n=\$(i+1);gsub(/"/, "", n);break}}
        print \$1,\$4-1,\$5,n,0,\$7
    }' ${gtf} | sort -k1,1 -k2,2n > gene.bed

    awk 'BEGIN{OFS="\\t"} !/^#/ && \$3=="exon" {
        n="."
        for(i=9;i<=NF;i++){if(\$i=="gene_id"){n=\$(i+1);gsub(/"/, "", n);break}}
        print \$1,\$4-1,\$5,n,0,\$7
    }' ${gtf} | sort -k1,1 -k2,2n > exon.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
