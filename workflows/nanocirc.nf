/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { NANOPLOT               } from '../modules/nf-core/nanoplot/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { FASTQ_RENAME           } from '../modules/local/fastq_rename'
include { CIRCRNA_ANALYSIS       } from '../subworkflows/local/circrna_analysis'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nanocirc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOCIRC {

    take:
    ch_samplesheet // channel: [ val(meta), path(fastq) ] from --input samplesheet

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // Ensure all FASTQ files are named *.fastq.gz — required by NanoPlot
    // which detects input type from extension
    FASTQ_RENAME ( ch_samplesheet )
    ch_fastq = FASTQ_RENAME.out.fastq

    //
    // QC: FastQC and NanoPlot
    //
    if (!params.skip_qc) {
        if (!params.skip_fastqc) {
            FASTQC ( ch_fastq )
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ it[1] })
            ch_versions      = ch_versions.mix(FASTQC.out.versions.first())
        }
        if (!params.skip_nanoplot) {
            NANOPLOT ( ch_fastq )
            ch_versions = ch_versions.mix(NANOPLOT.out.versions.first())
        }
    }

    //
    // SUBWORKFLOW: circRNA detection, BED12 conversion, merge and confidence scoring
    //
    if (!params.fasta) {
        error("CircRNA analysis requires '--fasta' (reference genome).")
    }
    if (!params.gtf) {
        error("CircRNA analysis requires '--gtf' (gene annotation).")
    }
    if ((params.run_isocirc || params.run_cirilong) && !params.circrna_db) {
        error("isocirc and ciri-long require '--circrna_db' (circRNA database).")
    }
    if (!params.run_isocirc && !params.run_circfl && !params.run_cirilong && !params.run_circnick) {
        error("At least one tool must be active. Use --run_isocirc, --run_circfl, --run_cirilong, or --run_circnick.")
    }
    if (params.run_circnick) {
        if (!params.circnick_species) {
            error("Parameter '--circnick_species' is required when '--run_circnick' is set. Valid options: 'mouse', 'human'")
        }
        if (params.circnick_species != 'mouse' && params.circnick_species != 'human') {
            error("Invalid --circnick_species: '${params.circnick_species}'. Valid options: 'mouse', 'human'")
        }
        if (!params.circnick_liftover_chain) {
            log.warn "No --circnick_liftover_chain provided. circnick-lrs will use built-in ${params.circnick_species == 'mouse' ? 'mm10' : 'hg19'} coordinates."
        }
    }

    CIRCRNA_ANALYSIS (
        ch_fastq,
        file(params.fasta, checkIfExists: true),
        file(params.gtf,   checkIfExists: true),
        params.circrna_db ? file(params.circrna_db, checkIfExists: true) : file('NO_FILE')
    )
    ch_versions = ch_versions.mix(CIRCRNA_ANALYSIS.out.versions)

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'nanocirc_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        ch_multiqc_config        = channel.fromPath(
            "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ?
            channel.fromPath(params.multiqc_config, checkIfExists: true) :
            channel.empty()
        ch_multiqc_logo          = params.multiqc_logo ?
            channel.fromPath(params.multiqc_logo, checkIfExists: true) :
            channel.empty()

        summary_params      = paramsSummaryMap(
            workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
            file(params.multiqc_methods_description, checkIfExists: true) :
            file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = channel.value(
            methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_methods_description.collectFile(
                name: 'methods_description_mqc.yaml',
                sort: true
            )
        )

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )
    }

    emit:
    multiqc_report = params.skip_multiqc ? [] : MULTIQC.out.report.toList()
    versions       = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
