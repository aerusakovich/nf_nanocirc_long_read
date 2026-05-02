/*
 * circRNA functional annotation subworkflow
 *
 * For each of the 6 merged outputs (strict/relaxed/exon × union/intersection):
 *   1. BED12  → GFF3          via agat_convert_bed2gff
 *   2. GFF3   → class codes   via gffcompare (against reference GTF)
 *   3. GFF3   → spliced FASTA via gffread    (against genome FASTA)
 *   4. TSV + .tmap → annotated TSV with class_code, ref_gene_id, ref_id appended
 *
 * Input channel shape:  [ val(meta), path(bed12), path(tsv) ]
 *   meta must contain a 'category' field (e.g. 'strict_union', 'exon_inter')
 *   used to build unique file prefixes and publish paths.
 *
 * Versions: all nf-core modules here use `topic: versions` and are collected
 * automatically by Channel.topic("versions") in the main workflow — no manual
 * ch_versions propagation is needed or correct for these modules.
 */

include { SAMTOOLS_FAIDX         } from '../../modules/nf-core/samtools/faidx/main'
include { AGAT_CONVERTBED2GFF    } from '../../modules/nf-core/agat/convertbed2gff/main'
include { GFFCOMPARE             } from '../../modules/nf-core/gffcompare/main'
include { GFFREAD                } from '../../modules/nf-core/gffread/main'
include { ADD_CLASS_CODES        } from '../../modules/local/add_class_codes'

workflow CIRCRNA_ANNOTATE {

    take:
    ch_beds_tsv  // channel: [ val(meta), path(bed12), path(tsv) ]
                 //          meta contains: id (sample), category (strict_union etc.)
    fasta        // path: genome FASTA
    gtf          // path: reference GTF

    main:

    // ── Index genome (once, reused for all samples / categories) ──────────
    SAMTOOLS_FAIDX(
        Channel.value([ [id: 'genome'], fasta, [] ]),
        false
    )

    // Build singleton genome channel [meta, fasta, fai] for GFFCOMPARE
    ch_genome = SAMTOOLS_FAIDX.out.fai
        .map { meta, fai -> [ meta, fasta, fai ] }

    // Singleton reference GTF channel for GFFCOMPARE
    ch_ref_gtf = Channel.value([ [id: 'ref'], gtf ])

    // ── Step 1: BED12 → GFF3 ─────────────────────────────────────────────
    AGAT_CONVERTBED2GFF(
        ch_beds_tsv.map { meta, bed, tsv -> [ meta, bed ] }
    )

    // ── Step 2: GFF3 + GTF → .tmap (class codes) ─────────────────────────
    GFFCOMPARE(
        AGAT_CONVERTBED2GFF.out.gff,
        ch_genome,
        ch_ref_gtf
    )

    // ── Step 3: GFF3 + genome FASTA → spliced sequences ──────────────────
    GFFREAD(
        AGAT_CONVERTBED2GFF.out.gff,
        fasta
    )

    // ── Step 4: Enrich confidence TSV with class codes ────────────────────
    // Join TSV + tmap + GFF so add_class_codes.py can resolve AGAT's integer
    // IDs back to BSJ IDs via the GFF Name attribute
    ch_tsv_tmap_gff = ch_beds_tsv
        .map  { meta, bed, tsv -> [ meta, tsv ] }
        .join ( GFFCOMPARE.out.tmap,           by: 0 )
        .join ( AGAT_CONVERTBED2GFF.out.gff,   by: 0 )

    ADD_CLASS_CODES(ch_tsv_tmap_gff)

    emit:
    annotated_gff = AGAT_CONVERTBED2GFF.out.gff   // [ meta, *.gff   ]
    spliced_fasta = GFFREAD.out.gffread_fasta       // [ meta, *.fasta ]
    annotated_tsv = ADD_CLASS_CODES.out.tsv         // [ meta, *.annotated.tsv ]
}
