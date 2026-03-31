/*
 * circRNA analysis subworkflow
 * Runs up to 4 detection tools, converts outputs to BED12,
 * and merges results when 2+ tools are active.
 */

include { ISOCIRC                } from '../../modules/local/isocirc'
include { CIRCFL_SEQ             } from '../../modules/local/circfl_seq'
include { PREPARE_GENOME         } from '../../modules/local/prepare_genome'
include { CIRI_LONG              } from '../../modules/local/ciri_long'
include { CIRCNICK_LRS           } from '../../modules/local/circnick_lrs'
include { CIRCNICK_LIFTOVER      } from '../../modules/local/circnick_liftover'
include { CIRILONG_TO_BED12      } from '../../modules/local/cirilong_to_bed12'
include { CIRCNICK_TO_BED12      } from '../../modules/local/circnick_to_bed12'
include { CIRCRNA_BEDTOOLS_PAIRS } from '../../modules/local/circrna_bedtools_pairs'
include { CIRCRNA_MERGE          } from '../../modules/local/circrna_merge'
include { CIRCRNA_EXON_MERGE     } from '../../modules/local/circrna_exon_merge'

workflow CIRCRNA_ANALYSIS {

    take:
    ch_fastq              // channel: [ val(meta), path(fastq) ]
    fasta                 // path: reference genome FASTA
    gtf                   // path: gene annotation GTF
    circrna_db            // path: circRNA database

    main:

    ch_versions = channel.empty()

    // ── STEP 1: run active tools ───────────────────────────────────────────

    // isocirc
    ch_isocirc_bed = channel.empty()
    if (params.run_isocirc) {
        ISOCIRC ( ch_fastq, fasta, gtf, circrna_db )
        ch_isocirc_bed = ISOCIRC.out.bed
        ch_versions    = ch_versions.mix(ISOCIRC.out.versions.first())
    }

    // circFL-seq
    ch_circfl_bed = channel.empty()
    if (params.run_circfl) {
        CIRCFL_SEQ ( ch_fastq, fasta, gtf )
        ch_circfl_bed = CIRCFL_SEQ.out.bed
        ch_versions   = ch_versions.mix(CIRCFL_SEQ.out.versions.first())
    }

    // CIRI-long
    ch_cirilong_bed = channel.empty()
    if (params.run_cirilong) {
        PREPARE_GENOME ( fasta )
        CIRI_LONG (
            ch_fastq,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.index,
            gtf,
            circrna_db
        )
        CIRILONG_TO_BED12 ( CIRI_LONG.out.info )
        ch_cirilong_bed = CIRILONG_TO_BED12.out.bed12
        ch_versions     = ch_versions.mix(PREPARE_GENOME.out.versions.first())
        ch_versions     = ch_versions.mix(CIRI_LONG.out.versions.first())
        ch_versions     = ch_versions.mix(CIRILONG_TO_BED12.out.versions.first())
    }

    // circnick-lrs
    ch_circnick_bed = channel.empty()
    if (params.run_circnick) {
        CIRCNICK_LRS ( ch_fastq, params.circnick_species )
        ch_versions = ch_versions.mix(CIRCNICK_LRS.out.versions.first())

        if (params.circnick_liftover_chain) {
            ch_chain = channel.fromPath(params.circnick_liftover_chain, checkIfExists: true)

            CIRCNICK_LIFTOVER (
                CIRCNICK_LRS.out.annotated,
                CIRCNICK_LRS.out.exon_usage,
                CIRCNICK_LRS.out.intron_cov,
                ch_chain.first()
            )
            ch_versions = ch_versions.mix(CIRCNICK_LIFTOVER.out.versions.first())

            CIRCNICK_TO_BED12 (
                CIRCNICK_LIFTOVER.out.annotated,
                CIRCNICK_LIFTOVER.out.exon_usage,
                CIRCNICK_LIFTOVER.out.intron_cov
            )
        } else {
            CIRCNICK_TO_BED12 (
                CIRCNICK_LRS.out.annotated,
                CIRCNICK_LRS.out.exon_usage,
                CIRCNICK_LRS.out.intron_cov
            )
        }

        ch_circnick_bed = CIRCNICK_TO_BED12.out.bed12
        ch_versions     = ch_versions.mix(CIRCNICK_TO_BED12.out.versions.first())
    }

    // ── STEP 2: collect active BED12s ─────────────────────────────────────
    ch_all_beds = channel.empty()

    if (params.run_isocirc) {
        ch_all_beds = ch_all_beds.mix(
            ch_isocirc_bed.map { meta, bed -> [ meta, 'isocirc', bed ] }
        )
    }
    if (params.run_circfl) {
        ch_all_beds = ch_all_beds.mix(
            ch_circfl_bed.map { meta, bed -> [ meta, 'circfl', bed ] }
        )
    }
    if (params.run_cirilong) {
        ch_all_beds = ch_all_beds.mix(
            ch_cirilong_bed.map { meta, bed -> [ meta, 'cirilong', bed ] }
        )
    }
    if (params.run_circnick) {
        ch_all_beds = ch_all_beds.mix(
            ch_circnick_bed.map { meta, bed -> [ meta, 'circnick', bed ] }
        )
    }

    // group by sample: [ meta, [tool_names], [bed_files] ]
    ch_beds_collected = ch_all_beds
        .groupTuple()
        .map { meta, tool_names, bed_files ->
            [ meta, tool_names, bed_files ]
        }

    // ── STEP 3: run merge pipeline only when 2+ tools active ──────────────
    ch_n_active = params.run_isocirc  ? 1 : 0
    ch_n_active = params.run_circfl   ? ch_n_active + 1 : ch_n_active
    ch_n_active = params.run_cirilong ? ch_n_active + 1 : ch_n_active
    ch_n_active = params.run_circnick ? ch_n_active + 1 : ch_n_active

    if (ch_n_active < 4) {
        log.warn(
            "[nf-core/nanocirc] Only ${ch_n_active} of 4 detection tools are active. " +
            "Tool consensus scores (Low/Medium/High) reflect agreement among the " +
            "tools that ran — a 'High' score from ${ch_n_active} tools is not " +
            "equivalent to 'High' from all 4 tools."
        )
    }

    if (ch_n_active >= 2) {

        CIRCRNA_BEDTOOLS_PAIRS ( ch_beds_collected )
        ch_versions = ch_versions.mix(CIRCRNA_BEDTOOLS_PAIRS.out.versions.first())

        ch_for_merge = ch_beds_collected
            .join(CIRCRNA_BEDTOOLS_PAIRS.out.pairs, by: 0)

        CIRCRNA_MERGE (
            ch_for_merge.map { meta, tool_names, bed_files, pairs -> [ meta, tool_names, bed_files ] },
            ch_for_merge.map { meta, tool_names, bed_files, pairs -> [ meta, pairs ] },
            ch_n_active
        )
        ch_versions = ch_versions.mix(CIRCRNA_MERGE.out.versions.first())

        CIRCRNA_EXON_MERGE (
            CIRCRNA_BEDTOOLS_PAIRS.out.pairs,
            ch_n_active
        )
        ch_versions = ch_versions.mix(CIRCRNA_EXON_MERGE.out.versions.first())
    }

    emit:
    isocirc_bed12   = ch_isocirc_bed
    circfl_bed12    = ch_circfl_bed
    cirilong_bed12  = ch_cirilong_bed
    circnick_bed12  = ch_circnick_bed

    strict_union_bed    = ch_n_active >= 2 ? CIRCRNA_MERGE.out.strict_union_bed    : channel.empty()
    strict_union_conf   = ch_n_active >= 2 ? CIRCRNA_MERGE.out.strict_union_conf   : channel.empty()
    strict_inter_bed    = ch_n_active >= 2 ? CIRCRNA_MERGE.out.strict_inter_bed    : channel.empty()
    strict_inter_conf   = ch_n_active >= 2 ? CIRCRNA_MERGE.out.strict_inter_conf   : channel.empty()
    relaxed_union_bed   = ch_n_active >= 2 ? CIRCRNA_MERGE.out.relaxed_union_bed   : channel.empty()
    relaxed_union_conf  = ch_n_active >= 2 ? CIRCRNA_MERGE.out.relaxed_union_conf  : channel.empty()
    relaxed_inter_bed   = ch_n_active >= 2 ? CIRCRNA_MERGE.out.relaxed_inter_bed   : channel.empty()
    relaxed_inter_conf  = ch_n_active >= 2 ? CIRCRNA_MERGE.out.relaxed_inter_conf  : channel.empty()
    exon_union_bed      = ch_n_active >= 2 ? CIRCRNA_EXON_MERGE.out.exon_union_bed  : channel.empty()
    exon_union_conf     = ch_n_active >= 2 ? CIRCRNA_EXON_MERGE.out.exon_union_conf : channel.empty()
    exon_inter_bed      = ch_n_active >= 2 ? CIRCRNA_EXON_MERGE.out.exon_inter_bed  : channel.empty()
    exon_inter_conf     = ch_n_active >= 2 ? CIRCRNA_EXON_MERGE.out.exon_inter_conf : channel.empty()

    versions = ch_versions
}
