process CIRI_LONG {
    tag "$meta.id"
    label 'process_high'

    containerOptions = "--writable-tmpfs"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/anrusakovich/ciri-long:latest' :
        'quay.io/anrusakovich/ciri-long:latest' }"

    input:
    tuple val(meta), path(fastq)
    path  fasta
    path  index
    path  gtf
    path  circrna_db

    output:
    tuple val(meta), path("${meta.id}_cirilong/${meta.id}.info"),       emit: info
    tuple val(meta), path("${meta.id}_cirilong/${meta.id}.isoforms"),   emit: isoforms
    tuple val(meta), path("${meta.id}_cirilong/"),                      emit: output_dir
    path  "versions.yml",                                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def out  = "${meta.id}_cirilong"
    """
    mkdir -p ${out}

    CIRI-long call \\
        -i ${fastq} \\
        -o ${out} \\
        -r ${fasta} \\
        -p ${meta.id} \\
        -a ${gtf} \\
        -c ${circrna_db} \\
        -t ${task.cpus} \\
        ${args}

    echo "${meta.id} ${out}/${meta.id}.cand_circ.fa" > ${out}/collapse_list.txt

    CIRI-long collapse \\
        -i ${out}/collapse_list.txt \\
        -o ${out} \\
        -r ${fasta} \\
        -p ${meta.id} \\
        -a ${gtf} \\
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciri_long: \$(CIRI-long --version 2>&1 | head -1)
    END_VERSIONS
    """
}
