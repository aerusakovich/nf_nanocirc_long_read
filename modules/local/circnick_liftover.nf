process CIRCNICK_LIFTOVER {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-liftover:377--ha8a8165_3' :
        'quay.io/biocontainers/ucsc-liftover:377--ha8a8165_3' }"

    input:
    tuple val(meta), path(annotated)    // circRNA_candidates.annotated.txt
    tuple val(meta), path(exon_usage)   // circ_circRNA_exon_usage_length_of_exons.txt
    tuple val(meta), path(intron_cov)   // intronCov.bed — also needs liftover
    path  chain                         // liftOver chain file (.chain or .chain.gz)

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
    # ── Step 1: lift annotated.txt ────────────────────────────
    tail -n +2 ${annotated} \\
        | awk 'OFS="\\t" {print \$2, \$3, \$4, \$1}' \\
        > annotated_for_lift.bed

    liftOver \\
        annotated_for_lift.bed \\
        ${chain} \\
        annotated_lifted.bed \\
        annotated_unmapped.bed

    # ── Step 2: lift exon_usage ───────────────────────────────
    tail -n +2 ${exon_usage} \\
        | awk 'OFS="\\t" {
            split(\$5, a, "_");
            chrom = a[length(a)];
            split(chrom, b, ":");
            print b[1], \$6, \$7, \$1 "_EXON_" NR
        }' \\
        > exon_for_lift.bed

    liftOver \\
        exon_for_lift.bed \\
        ${chain} \\
        exon_lifted.bed \\
        exon_unmapped.bed

    # ── Step 3: lift intron_cov.bed ───────────────────────────
    # intron_cov is already BED format — cols 0,1,2 are chrom,start,end
    # col 10 is circ_id — keep as name for re-joining
    awk 'OFS="\\t" {print \$1, \$2, \$3, \$11 "_INTRON_" NR}' ${intron_cov} \\
        > intron_for_lift.bed

    liftOver \\
        intron_for_lift.bed \\
        ${chain} \\
        intron_lifted.bed \\
        intron_unmapped.bed

    # ── Step 4: rejoin all lifted coords ─────────────────────
    python3 ${projectDir}/bin/circnick_liftover.py \\
        --lifted_annotated  annotated_lifted.bed \\
        --lifted_exons      exon_lifted.bed \\
        --orig_annotated    ${annotated} \\
        --orig_exon_usage   ${exon_usage} \\
        --sample            ${meta.id}

    # ── Step 5: rebuild lifted intron_cov from lifted coords ─
    python3 ${projectDir}/bin/circnick_liftover_introns.py \\
        --lifted_introns  intron_lifted.bed \\
        --orig_introns    ${intron_cov} \\
        --sample          ${meta.id}

    echo "Unmapped annotated : \$(grep -v '^#' annotated_unmapped.bed | wc -l)"
    echo "Unmapped exons     : \$(grep -v '^#' exon_unmapped.bed | wc -l)"
    echo "Unmapped introns   : \$(grep -v '^#' intron_unmapped.bed | wc -l)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        liftover: \$(liftOver 2>&1 | head -1 | sed 's/liftOver - //')
    END_VERSIONS
    """
}
