process CIRCNICK_LRS {
    tag "$meta.id"
    label 'process_high'

    containerOptions = "--writable-tmpfs"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/anrusakovich/circnick-lrs:latest' :
        'quay.io/anrusakovich/circnick-lrs:latest' }"

    input:
    tuple val(meta), path(fastq)
    val   species

    output:
    tuple val(meta), path("${meta.id}_circnick/${meta.id}/${meta.id}.circRNA_candidates.annotated.txt"),                                          emit: annotated
    tuple val(meta), path("${meta.id}_circnick/${meta.id}/${meta.id}.circ_circRNA_exon_usage_length_of_exons.txt"),                               emit: exon_usage
    tuple val(meta), path("${meta.id}_circnick/${meta.id}/${meta.id}.introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.intronCov.bed"),     emit: intron_cov
    tuple val(meta), path("${meta.id}_circnick/"),                                                                                                emit: output_dir
    path  "versions.yml",                                                                                                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def out   = "${meta.id}_circnick"
    def is_gz = fastq.name.endsWith('.gz')
    """
    mkdir -p ${out}

    # circnick-lrs requires gzipped input.
    # Use meta.id as filename — circnick uses it as the output subdirectory name.
    ${is_gz ?
        "ln -sf ${fastq} ${meta.id}.fq.gz" :
        "gzip -c ${fastq} > ${meta.id}.fq.gz"
    }

    long_read_circRNA run \\
        --species          ${species} \\
        --reference-path   /opt/long_read_circRNA/data \\
        --script-path      /opt/long_read_circRNA/scripts \\
        --output-path      ${out} \\
        ${args} \\
        ${meta.id}.fq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circnick_lrs: \$(grep '__version__' /opt/long_read_circRNA/long_read_circRNA | head -1 | sed 's/.*= *//')
    END_VERSIONS
    """
}
