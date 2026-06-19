process CIRCNICK_LIFTOVER {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/liftover:1.4.0--py313h7fbb527_0' :
        'quay.io/biocontainers/liftover:1.4.0--py313h7fbb527_0' }"

    input:
    tuple val(meta),  path(annotated)    // circRNA_candidates.annotated.txt
    tuple val(meta2), path(exon_usage)   // circ_circRNA_exon_usage_length_of_exons.txt
    tuple val(meta3), path(intron_cov)   // intronCov.bed — also needs liftover
    path  chain                          // liftOver chain file (.chain or .chain.gz)

    output:
    tuple val(meta), path("${meta.id}_lifted_annotated.txt"),  emit: annotated
    tuple val(meta), path("${meta.id}_lifted_exon_usage.txt"), emit: exon_usage
    tuple val(meta), path("${meta.id}_lifted_intron_cov.bed"), emit: intron_cov
    tuple val(meta), path("${meta.id}_liftover_failed.tsv"),   emit: failed
    path  "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export HOME=\$PWD

    # ── Step 1: lift annotated.txt ────────────────────────────
    tail -n +2 ${annotated} \\
        | awk 'OFS="\\t" {print \$2, \$3, \$4, \$1}' \\
        > annotated_for_lift.bed

    liftover_bed.py \\
        annotated_for_lift.bed \\
        ${chain} \\
        annotated_lifted.bed \\
        annotated_unmapped.bed

    # ── Step 2: lift exon_usage ───────────────────────────────
    tail -n +2 ${exon_usage} \\
        | awk 'BEGIN{OFS="\\t"} {
            split(\$5, a, "_"); chrom = a[length(a)];
            if (index(chrom, ":") > 0) {
                split(chrom, b, ":");
                print b[1], \$6, \$7, \$1 "_EXON_" NR
            }
        }' \\
        > exon_for_lift.bed

    liftover_bed.py \\
        exon_for_lift.bed \\
        ${chain} \\
        exon_lifted.bed \\
        exon_unmapped.bed

    # ── Step 3: lift intron_cov.bed ───────────────────────────
    # intron_cov is already BED format — cols 0,1,2 are chrom,start,end
    # col 10 is circ_id — keep as name for re-joining
    awk 'OFS="\\t" {print \$1, \$2, \$3, \$11 "_INTRON_" NR}' ${intron_cov} \\
        > intron_for_lift.bed

    liftover_bed.py \\
        intron_for_lift.bed \\
        ${chain} \\
        intron_lifted.bed \\
        intron_unmapped.bed

    # ── Step 4: rejoin all lifted coords ─────────────────────
    circnick_liftover.py \\
        --lifted_annotated  annotated_lifted.bed \\
        --lifted_exons      exon_lifted.bed \\
        --orig_annotated    ${annotated} \\
        --orig_exon_usage   ${exon_usage} \\
        --sample            ${meta.id}

    # ── Step 5: rebuild lifted intron_cov from lifted coords ─
    circnick_liftover_introns.py \\
        --lifted_introns  intron_lifted.bed \\
        --orig_introns    ${intron_cov} \\
        --sample          ${meta.id}

    echo "Unmapped annotated : \$(grep -v '^#' annotated_unmapped.bed | wc -l)"
    echo "Unmapped exons     : \$(grep -v '^#' exon_unmapped.bed | wc -l)"
    echo "Unmapped introns   : \$(grep -v '^#' intron_unmapped.bed | wc -l)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        liftover: \$(python -c "import liftover; print(liftover.__version__)")
    END_VERSIONS
    """
}
