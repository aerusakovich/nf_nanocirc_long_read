process CIRCFL_SEQ {
    tag "$meta.id"
    label 'process_high'

    containerOptions = "--writable-tmpfs"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/anrusakovich/circfl-seq:latest' :
        'quay.io/anrusakovich/circfl-seq:latest' }"

    input:
    tuple val(meta), path(fastq)
    path  fasta
    path  gtf

    output:
    tuple val(meta), path("${meta.id}_circfl/circFL_final.bed"), emit: bed
    tuple val(meta), path("${meta.id}_circfl/"),                 emit: output_dir
    path  "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def out   = "${meta.id}_circfl"
    def is_gz = fastq.name.endsWith('.gz')
    """
    export PYTHONNOUSERSITE=1

    # Decompress inline if needed — circfull requires plain FASTQ
    ${is_gz ? "gunzip -c ${fastq} > input.fastq" : "ln -sf ${fastq} input.fastq"}

    # Sort and index GTF
    grep "^#" ${gtf} > anno_sorted.gtf
    grep -v "^#" ${gtf} | sort -k1,1 -k4,4n >> anno_sorted.gtf
    bgzip -c anno_sorted.gtf > anno.gtf.gz
    tabix -p gff anno.gtf.gz

    mkdir -p ${out}

    circfull RG \\
        -f input.fastq \\
        -g ${fasta} \\
        -a anno.gtf.gz \\
        -t ${task.cpus} \\
        -r \\
        -o ${out}

    circfull DNSC \\
        -f ${out}/RG/circSeq.fa \\
        -t ${task.cpus} \\
        -o ${out}

    circfull cRG \\
        -f ${out}/DNSC \\
        -g ${fasta} \\
        -a anno.gtf.gz \\
        -t ${task.cpus} \\
        -o ${out}

    circfull mRG \\
        -f input.fastq \\
        -g ${fasta} \\
        -r ${out} \\
        -c ${out} \\
        -o ${out}

    circfull anno \\
        -b ${out}/mRG/circFL_Normal_pass.bed \\
        -a anno.gtf.gz \\
        -o ${out}/annotated

    circfull FL2BED \\
        -c ${out}/mRG/circFL_Normal_pass.txt \\
        -o ${out}/circFL_final.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circfl_seq: \$(circfull --version 2>&1 | head -1)
    END_VERSIONS
    """
}
