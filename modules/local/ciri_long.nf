process CIRI_LONG {
    tag "$meta.id"
    label 'process_high'

    containerOptions "--writable-tmpfs"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://zenodo.org/records/20707975/files/nanocirc-ciri-long-v1.0.sif?download=1' :
        'quay.io/anrusakovich/ciri-long:latest' }"

    input:
    tuple val(meta), path(fastq)
    path  fasta
    path  index
    path  gtf
    path  circrna_db

    output:
    tuple val(meta), path("${meta.id}_cirilong/${meta.id}.info"),         emit: info
    tuple val(meta), path("${meta.id}_cirilong/${meta.id}.isoforms"),     emit: isoforms
    tuple val(meta), path("${meta.id}_cirilong/${meta.id}.expression", optional: true),   emit: expr
    tuple val(meta), path("${meta.id}_cirilong/"),                        emit: output_dir
    path  "versions.yml",                                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def out  = "${meta.id}_cirilong"
    """
    mkdir -p ${out}

    # CIRI-long requires the -c file to:
    #   (a) have a .bed suffix — it checks circ_path.suffix == ".bed"
    #   (b) be BED4: chrom  start  end  strand  (col[3] is read as strand)
    # Sanitise whatever format was passed: strip to BED4 and force .bed extension.
    python3 - <<'PYEOF'
    import sys
    out_path = "circrna_db_clean.bed"
    with open("${circrna_db}") as fin, open(out_path, "w") as fout:
        for line in fin:
            c = line.rstrip("\\n").split("\\t")
            if len(c) < 4:
                continue
            # BED4 → strand already in col[3]
            # BED6+ → strand is in col[5]
            strand = c[5] if len(c) >= 6 else c[3]
            fout.write("\\t".join([c[0], c[1], c[2], strand]) + "\\n")
    PYEOF

    # Remove any truncated tmp files left by a previous crashed run.
    # CIRI-long reuses tmp/{prefix}.ccs.fa / tmp/{prefix}.raw.fa when they exist;
    # if those files are out of sync (crash mid-write) the Python worker raises
    # "ValueError: not enough values to unpack (expected 4, got 3)".
    if [ -d "${out}/tmp" ]; then
        ccs="${out}/tmp/${meta.id}.ccs.fa"
        raw="${out}/tmp/${meta.id}.raw.fa"
        if [ -f "\$ccs" ] && [ -f "\$raw" ]; then
            n_ccs=\$(grep -c '^>' "\$ccs" || true)
            n_raw=\$(grep -c '^>' "\$raw" || true)
            if [ "\$n_ccs" != "\$n_raw" ]; then
                echo "[ciri_long] tmp mismatch (ccs=\$n_ccs raw=\$n_raw) — removing to force fresh extraction"
                rm -f "\$ccs" "\$raw"
            fi
        fi
    fi

    CIRI-long call \\
        -i ${fastq} \\
        -o ${out} \\
        -r ${fasta} \\
        -p ${meta.id} \\
        -a ${gtf} \\
        -c circrna_db_clean.bed \\
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
