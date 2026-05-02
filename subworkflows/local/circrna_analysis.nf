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
include { CIRCRNA_SMART_MERGE    } from '../../modules/local/circrna_smart_merge'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_BALANCED           } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_HIGH_CONFIDENCE    } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_CONSENSUS_NO_LOW   } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_CONSENSUS_TRUSTED  } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_XSTRUCT_NO_LOW     } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_XSTRUCT_TRUSTED    } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_PRIORITY_NO_LOW    } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_PRIORITY_TRUSTED   } from '../../modules/local/circrna_confidence_filter'
include { GTF_TO_FEATURE_BED     } from '../../modules/local/gtf_to_feature_bed'
include { CIRCRNA_FINALIZE       } from '../../modules/local/circrna_finalize'
include { CIRCRNA_ANNOTATE       } from './circrna_annotate'

workflow CIRCRNA_ANALYSIS {

    take:
    ch_fastq              // channel: [ val(meta), path(fastq) ]
    fasta                 // path: reference genome FASTA
    gtf                   // path: gene annotation GTF
    circrna_db            // path: circRNA database

    main:

    ch_versions = channel.empty()

    // Gene/exon BED files derived from GTF — used for type classification
    GTF_TO_FEATURE_BED(gtf)
    ch_versions = ch_versions.mix(GTF_TO_FEATURE_BED.out.versions)
    ch_gene_bed = GTF_TO_FEATURE_BED.out.gene_bed.first()
    ch_exon_bed = GTF_TO_FEATURE_BED.out.exon_bed.first()

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

    // Group by sample: [ meta, [tool_names], [bed_files] ]
    ch_beds_collected = ch_all_beds
        .groupTuple()
        .map { meta, tool_names, bed_files ->
            [ meta, tool_names, bed_files ]
        }

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

        // Legacy merge modes (--run_benchmark_modes)
        if (params.run_benchmark_modes) {
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

        CIRCRNA_SMART_MERGE (
            ch_for_merge.map { meta, tool_names, bed_files, pairs -> [ meta, tool_names, bed_files ] },
            ch_for_merge.map { meta, tool_names, bed_files, pairs -> [ meta, pairs ] },
            ch_n_active
        )
        ch_versions = ch_versions.mix(CIRCRNA_SMART_MERGE.out.versions.first())

        // Build per-mode filter input channels — both use hybrid as input
        ch_hybrid = CIRCRNA_SMART_MERGE.out.hybrid_bed
            .join(CIRCRNA_SMART_MERGE.out.hybrid_conf, by: 0)

        ch_for_balanced        = ch_hybrid.map { meta, bed, tsv -> [ meta + [category: 'hybrid'], bed, tsv ] }
        ch_for_high_confidence = ch_hybrid.map { meta, bed, tsv -> [ meta + [category: 'hybrid'], bed, tsv ] }

        CIRCRNA_FILTER_BALANCED        ( ch_for_balanced )
        CIRCRNA_FILTER_HIGH_CONFIDENCE ( ch_for_high_confidence )
        ch_versions = ch_versions.mix(CIRCRNA_FILTER_BALANCED.out.versions.first())

        if (params.run_benchmark_modes) {
            ch_consensus = CIRCRNA_SMART_MERGE.out.consensus_bed
                .join(CIRCRNA_SMART_MERGE.out.consensus_conf, by: 0)

            ch_xstruct = CIRCRNA_SMART_MERGE.out.consensus_xstruct_bed
                .join(CIRCRNA_SMART_MERGE.out.consensus_xstruct_conf, by: 0)

            ch_for_priority = CIRCRNA_SMART_MERGE.out.priority_bed
                .join(CIRCRNA_SMART_MERGE.out.priority_conf, by: 0)

            CIRCRNA_FILTER_CONSENSUS_NO_LOW  ( ch_consensus.map { meta, bed, tsv -> [ meta + [category: 'smart_consensus_no_low'],      bed, tsv ] } )
            CIRCRNA_FILTER_CONSENSUS_TRUSTED ( ch_consensus.map { meta, bed, tsv -> [ meta + [category: 'smart_consensus_filtered'],    bed, tsv ] } )
            CIRCRNA_FILTER_XSTRUCT_NO_LOW    ( ch_xstruct.map   { meta, bed, tsv -> [ meta + [category: 'smart_consensus_xstruct_no_low'],  bed, tsv ] } )
            CIRCRNA_FILTER_XSTRUCT_TRUSTED   ( ch_xstruct.map   { meta, bed, tsv -> [ meta + [category: 'smart_consensus_xstruct_filtered'], bed, tsv ] } )
            CIRCRNA_FILTER_PRIORITY_NO_LOW   ( ch_for_priority.map { meta, bed, tsv -> [ meta + [category: 'smart_priority_no_low'],    bed, tsv ] } )
            CIRCRNA_FILTER_PRIORITY_TRUSTED  ( ch_for_priority.map { meta, bed, tsv -> [ meta + [category: 'smart_priority_filtered'],  bed, tsv ] } )
        }

        // Discovery is hybrid emitted directly (no filter applied)
        ch_discovery_for_annotate = ch_hybrid
            .map { meta, bed, tsv -> [ meta + [category: 'discovery'], bed, tsv ] }

        if (!params.skip_annotation) {

            ch_for_annotate = ch_discovery_for_annotate
                .mix( CIRCRNA_FILTER_BALANCED.out.bed       .join(CIRCRNA_FILTER_BALANCED.out.conf,        by: 0).map { meta, bed, tsv -> [ meta + [category: 'balanced'],        bed, tsv ] } )
                .mix( CIRCRNA_FILTER_HIGH_CONFIDENCE.out.bed.join(CIRCRNA_FILTER_HIGH_CONFIDENCE.out.conf, by: 0).map { meta, bed, tsv -> [ meta + [category: 'high_confidence'], bed, tsv ] } )

            if (params.run_benchmark_modes) {
                ch_for_annotate = ch_for_annotate
                    .mix( CIRCRNA_MERGE.out.strict_union_bed  .join(CIRCRNA_MERGE.out.strict_union_conf,  by: 0).map { meta, bed, tsv -> [ meta + [category: 'strict_union'],  bed, tsv ] } )
                    .mix( CIRCRNA_MERGE.out.strict_inter_bed  .join(CIRCRNA_MERGE.out.strict_inter_conf,  by: 0).map { meta, bed, tsv -> [ meta + [category: 'strict_inter'],  bed, tsv ] } )
                    .mix( CIRCRNA_MERGE.out.relaxed_union_bed .join(CIRCRNA_MERGE.out.relaxed_union_conf, by: 0).map { meta, bed, tsv -> [ meta + [category: 'relaxed_union'], bed, tsv ] } )
                    .mix( CIRCRNA_MERGE.out.relaxed_inter_bed .join(CIRCRNA_MERGE.out.relaxed_inter_conf, by: 0).map { meta, bed, tsv -> [ meta + [category: 'relaxed_inter'], bed, tsv ] } )
                    .mix( CIRCRNA_EXON_MERGE.out.exon_union_bed.join(CIRCRNA_EXON_MERGE.out.exon_union_conf, by: 0).map { meta, bed, tsv -> [ meta + [category: 'exon_union'], bed, tsv ] } )
                    .mix( CIRCRNA_EXON_MERGE.out.exon_inter_bed.join(CIRCRNA_EXON_MERGE.out.exon_inter_conf, by: 0).map { meta, bed, tsv -> [ meta + [category: 'exon_inter'], bed, tsv ] } )
            }

            CIRCRNA_ANNOTATE(ch_for_annotate, fasta, gtf)
            // CIRCRNA_ANNOTATE uses nf-core modules with `topic: versions` —
            // versions are collected automatically; no ch_versions mixing needed.

            // ── Finalize: type + expression + clean TSV ─────────────────────────────
            // annotated_tsv has N items per sample (one per category); expression files have 1.
            // Use combine(by:0) for N:1 matching, not join() which is 1:1 only.
            // circfl emits a list (mRG + RG pass files); pick mRG when present.
            // ciri_long uses optional:true — fill missing samples with NO_FILE via remainder join.
            ch_iso_expr  = params.run_isocirc  ? ISOCIRC.out.expr.map { m, f -> [m.id, f] }
                                               : ch_fastq.map { m, _ -> [m.id, file('NO_FILE')] }
            ch_fl_expr   = params.run_circfl   ? CIRCFL_SEQ.out.expr.map { m, files ->
                               def fl = files instanceof List ? files : [files]
                               [m.id, fl.find { it.toString().contains('/mRG/') } ?: fl[0]]
                           }                   : ch_fastq.map { m, _ -> [m.id, file('NO_FILE')] }
            ch_nick_expr = params.run_circnick ? CIRCNICK_LRS.out.annotated.map { m, f -> [m.id, f] }
                                               : ch_fastq.map { m, _ -> [m.id, file('NO_FILE')] }
            // ciri_long: guarantee one entry per sample even when expression file is absent
            ch_ciri_expr = params.run_cirilong
                ? ch_fastq.map { m, _ -> [m.id] }
                    .join(CIRI_LONG.out.expr.map { m, f -> [m.id, f] }, remainder: true)
                    .map { id, f -> [id, f ?: file('NO_FILE')] }
                : ch_fastq.map { m, _ -> [m.id, file('NO_FILE')] }

            ch_for_finalize = CIRCRNA_ANNOTATE.out.annotated_tsv
                .map     { meta, tsv -> [meta.id, meta, tsv] }
                .combine ( ch_iso_expr,  by: 0 )
                .combine ( ch_fl_expr,   by: 0 )
                .combine ( ch_nick_expr, by: 0 )
                .combine ( ch_ciri_expr, by: 0 )
                .map     { id, meta, tsv, iso, fl, nick, ciri -> [meta, tsv, iso, ciri, fl, nick] }

            CIRCRNA_FINALIZE(ch_for_finalize, ch_gene_bed, ch_exon_bed)
            ch_versions = ch_versions.mix(CIRCRNA_FINALIZE.out.versions.first())
        }
    }

    emit:
    // Per-tool BED12 outputs (always present when tool is active)
    isocirc_bed12   = ch_isocirc_bed
    circfl_bed12    = ch_circfl_bed
    cirilong_bed12  = ch_cirilong_bed
    circnick_bed12  = ch_circnick_bed

    // Feature BED files (for type classification, passed to crossrun merge)
    gene_bed = ch_gene_bed
    exon_bed = ch_exon_bed

    // Discovery: hybrid, all isoforms, unfiltered
    discovery_bed  = ch_n_active >= 2 ? CIRCRNA_SMART_MERGE.out.hybrid_bed  : channel.empty()
    discovery_conf = ch_n_active >= 2 ? CIRCRNA_SMART_MERGE.out.hybrid_conf : channel.empty()
    // Balanced: hybrid + no_low filter
    balanced_bed   = ch_n_active >= 2 ? CIRCRNA_FILTER_BALANCED.out.bed        : channel.empty()
    balanced_conf  = ch_n_active >= 2 ? CIRCRNA_FILTER_BALANCED.out.conf       : channel.empty()
    // High-confidence: hybrid + high_only filter
    high_conf_bed  = ch_n_active >= 2 ? CIRCRNA_FILTER_HIGH_CONFIDENCE.out.bed  : channel.empty()
    high_conf_conf = ch_n_active >= 2 ? CIRCRNA_FILTER_HIGH_CONFIDENCE.out.conf : channel.empty()

    // Legacy merge modes (--run_benchmark_modes)
    strict_union_bed    = (ch_n_active >= 2 && params.run_benchmark_modes) ? CIRCRNA_MERGE.out.strict_union_bed    : channel.empty()
    strict_union_conf   = (ch_n_active >= 2 && params.run_benchmark_modes) ? CIRCRNA_MERGE.out.strict_union_conf   : channel.empty()
    strict_inter_bed    = (ch_n_active >= 2 && params.run_benchmark_modes) ? CIRCRNA_MERGE.out.strict_inter_bed    : channel.empty()
    strict_inter_conf   = (ch_n_active >= 2 && params.run_benchmark_modes) ? CIRCRNA_MERGE.out.strict_inter_conf   : channel.empty()
    relaxed_union_bed   = (ch_n_active >= 2 && params.run_benchmark_modes) ? CIRCRNA_MERGE.out.relaxed_union_bed   : channel.empty()
    relaxed_union_conf  = (ch_n_active >= 2 && params.run_benchmark_modes) ? CIRCRNA_MERGE.out.relaxed_union_conf  : channel.empty()
    relaxed_inter_bed   = (ch_n_active >= 2 && params.run_benchmark_modes) ? CIRCRNA_MERGE.out.relaxed_inter_bed   : channel.empty()
    relaxed_inter_conf  = (ch_n_active >= 2 && params.run_benchmark_modes) ? CIRCRNA_MERGE.out.relaxed_inter_conf  : channel.empty()
    exon_union_bed      = (ch_n_active >= 2 && params.run_benchmark_modes) ? CIRCRNA_EXON_MERGE.out.exon_union_bed  : channel.empty()
    exon_union_conf     = (ch_n_active >= 2 && params.run_benchmark_modes) ? CIRCRNA_EXON_MERGE.out.exon_union_conf : channel.empty()
    exon_inter_bed      = (ch_n_active >= 2 && params.run_benchmark_modes) ? CIRCRNA_EXON_MERGE.out.exon_inter_bed  : channel.empty()
    exon_inter_conf     = (ch_n_active >= 2 && params.run_benchmark_modes) ? CIRCRNA_EXON_MERGE.out.exon_inter_conf : channel.empty()

    annotated_gff = (!params.skip_annotation && ch_n_active >= 2) ? CIRCRNA_ANNOTATE.out.annotated_gff : channel.empty()
    spliced_fasta = (!params.skip_annotation && ch_n_active >= 2) ? CIRCRNA_ANNOTATE.out.spliced_fasta : channel.empty()
    annotated_tsv = (!params.skip_annotation && ch_n_active >= 2) ? CIRCRNA_ANNOTATE.out.annotated_tsv : channel.empty()
    clean_tsv     = (!params.skip_annotation && ch_n_active >= 2) ? CIRCRNA_FINALIZE.out.clean          : channel.empty()

    versions = ch_versions
}
