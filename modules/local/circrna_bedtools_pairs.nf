process CIRCRNA_BEDTOOLS_PAIRS {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), val(tool_names), path(bed_files)
    // tool_names : list of tool name strings e.g. ['isocirc', 'circfl', 'cirilong']
    // bed_files  : list of BED12 paths in the same order as tool_names

    output:
    tuple val(meta), path("*_vs_*.txt"), emit: pairs
    path  "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def names_str = tool_names.join(" ")
    def files_str = bed_files.collect { it.toString() }.join(" ")
    """
    tool_names=(${names_str})
    bed_files=(${files_str})

    n=\${#tool_names[@]}

    for (( i=0; i<n; i++ )); do
        for (( j=i+1; j<n; j++ )); do
            name_a=\${tool_names[i]}
            name_b=\${tool_names[j]}
            file_a=\${bed_files[i]}
            file_b=\${bed_files[j]}

            echo "Running bedtools: \${name_a} vs \${name_b}"

            bedtools intersect \\
                -a \${file_a} \\
                -b \${file_b} \\
                -split \\
                -wo \\
                > \${name_a}_vs_\${name_b}.txt

            echo "  \$(wc -l < \${name_a}_vs_\${name_b}.txt) overlapping pairs found"
        done
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools //')
    END_VERSIONS
    """
}
