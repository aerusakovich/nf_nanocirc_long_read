process PREPARE_GENOME {
    tag "prepare_genome"
    label 'process_low'

    // Cache results — skipped on subsequent runs if files already exist here
    storeDir "${params.genome_index_dir ?: "${params.outdir}/genome_index"}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwakit:0.7.15--he513fc3_2' :
        'quay.io/biocontainers/bwakit:0.7.15--he513fc3_2' }"

    input:
    path fasta

    output:
    path "genome.fa",    emit: fasta
    path "genome.fa.*",  emit: index
    path "versions.yml", emit: versions

    script:
    """
    # Resolve symlink — bwapy needs a real file to find index alongside it
    cp -L ${fasta} genome.fa

    # Copy index if it already exists next to source FASTA, otherwise build
    if [ -f "${fasta}.amb" ]; then
        echo "Index found next to source FASTA — copying"
        cp -L ${fasta}.amb genome.fa.amb
        cp -L ${fasta}.ann genome.fa.ann
        cp -L ${fasta}.bwt genome.fa.bwt
        cp -L ${fasta}.pac genome.fa.pac
        cp -L ${fasta}.sa  genome.fa.sa
    else
        echo "Index not found — building with bwa index"
        bwa index genome.fa
    fi

    if [ -f "${fasta}.fai" ]; then
        cp -L ${fasta}.fai genome.fa.fai
    else
        samtools faidx genome.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep Version | sed 's/Version: //')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
