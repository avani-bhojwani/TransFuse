/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowTransfuse.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                 } from '../subworkflows/local/input_check'
include { TRFORMAT                    } from '../modules/local/trformat'
include { TR2AACDS                    } from '../modules/local/tr2aacds'
include { RNAQUAST                    } from '../modules/local/rnaquast'
include { RNASPADES                   } from '../modules/local/rnaspades'
include { ASSEMBLE                    } from '../subworkflows/local/assemble'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTP                       } from '../modules/nf-core/fastp/main'
include { STAR_GENOMEGENERATE         } from '../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN                  } from '../modules/nf-core/star/align/main'
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { TRINITY                     } from '../modules/nf-core/trinity/main'
include { SALMON_INDEX                } from '../modules/nf-core/salmon/index/main'
include { SALMON_QUANT                } from '../modules/nf-core/salmon/quant/main'
include { BUSCO; BUSCO as BUSCO_COMBINED } from '../modules/nf-core/busco/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow TRANSFUSE {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // MODULE: FASTP
    //
    FASTP (
        INPUT_CHECK.out.reads,
        params.adapter_fasta,
        params.save_trimmed_fail,
        params.save_merged
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    // Method 1 pools all reads together and assembles them
    if (params.assemble_method == '1') {
        method_1_pool_ch = FASTP.out.reads.collect { meta, fastq -> fastq }.map {[ [id:'all_samples', single_end:false], it ] }

        //
        // MODULE: CAT_FASTQ
        //
        CAT_FASTQ (
            method_1_pool_ch
        )

        //
        // MODULE: ASSEMBLE
        //
        ASSEMBLE (
            CAT_FASTQ.out.reads
        )
        ch_versions = ch_versions.mix(ASSEMBLE.out.versions)
    }

    //when mapping to reference transcriptome
    if (!params.skip_mapping) {
        // 
        // MODULE: STAR genomeGenerate
        //
        STAR_GENOMEGENERATE (
            params.fasta,
            params.gtf,
            params.star_ignore_sjdbgtf
        )
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

        //
        // MODULE: STAR align
        //
        STAR_ALIGN (
            FASTP.out.reads,
            STAR_GENOMEGENERATE.out.index,
            params.gtf,
            params.star_ignore_sjdbgtf,
            params.seq_platform,
            params.seq_center
        )
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

        //
        // MODULE: Trinity
        //
        TRINITY (
            STAR_ALIGN.out.fastq
        )
        ch_versions = ch_versions.mix(TRINITY.out.versions.first())

        //
        // MODULE: RNA-SPAdes
        //
        RNASPADES (
            params.kmers,
            STAR_ALIGN.out.fastq
        )
        ch_version = ch_versions.mix(RNASPADES.out.versions)

        //
        // MODULE TRFORMAT
        //
        TRFORMAT (
            params.fasta,
            TRINITY.out.trinity_assembly,
            RNASPADES.out.spades_assembly
        )

        //
        // MODULE: BUSCO (for old reference transcriptome)
        //
        old_ref_ch = Channel.of(['old', params.fasta])
        BUSCO (
            old_ref_ch,
            params.busco_lineage,
            params.busco_lineages_path,
            params.busco_config
        )
        ch_versions = ch_versions.mix(BUSCO.out.versions)
    }

    //when skipping mapping to reference transcriptome
    else {
        //
        // MODULE: Trinity
        //
        TRINITY (
            FASTP.out.reads
        )
        ch_versions = ch_versions.mix(TRINITY.out.versions.first())

        //
        // MODULE: RNA-SPAdes
        //
        RNASPADES (
            params.kmers,
            FASTP.out.reads
        )
        ch_version = ch_versions.mix(RNASPADES.out.versions)

        //
        // MODULE TRFORMAT
        //
        TRFORMAT (
            [],
            TRINITY.out.trinity_assembly,
            RNASPADES.out.spades_assembly
        )
    }

    //
    // MODULE TR2AACDS
    //
    TR2AACDS (
        TRFORMAT.out.reformatted_fasta
    )
    ch_versions = ch_versions.mix(TR2AACDS.out.versions)

    //
    // MODULE: BUSCO combined (for newly assembled reference transcriptome)
    //
    new_ref_ch = TR2AACDS.out.non_redundant_fasta.map{ file -> ['new', file]}
    BUSCO_COMBINED (
        new_ref_ch,
        params.busco_lineage,
        params.busco_lineages_path,
        params.busco_config
    )

    //
    // MODULE: RNAQUAST
    //
    RNAQUAST (
        TR2AACDS.out.non_redundant_fasta
    )
    ch_versions = ch_versions.mix(RNAQUAST.out.versions)

    //
    // MODULE: Salmon index
    //
    SALMON_INDEX (
        TR2AACDS.out.non_redundant_fasta
    )
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)

    //
    // MODULE: Salmon quant
    //
    SALMON_QUANT (
        FASTP.out.reads,
        SALMON_INDEX.out.index,
        TR2AACDS.out.non_redundant_fasta,
        params.lib_type        
    )
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary = WorkflowTransfuse.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description = WorkflowTransfuse.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))

    if (!params.skip_mapping) {
        ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(BUSCO.out.short_summaries_txt.collect{it[1]}.ifEmpty([]))
    }

    ch_multiqc_files = ch_multiqc_files.mix(BUSCO_COMBINED.out.short_summaries_txt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.results.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
