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

// Check rRNA databases for sortmerna
if (params.remove_ribo_rna) {
    ch_ribo_db = file(params.ribo_database_manifest, checkIfExists: true)
    if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
}

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
// MODULE: Loaded from modules/local/
//
include { CAT_FASTA                   } from '../modules/local/cat_fasta'
include { MERGE_TPM                   } from '../modules/local/merge_tpm'
include { RNAQUAST                    } from '../modules/local/rnaquast'
include { RNASPADES                   } from '../modules/local/rnaspades'
include { TR2AACDS                    } from '../modules/local/tr2aacds'

//
// SUBWORKFLOW: Loaded from subworkflows/local/
//
include { INPUT_CHECK                 } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core
//
include { BUSCO                       } from '../modules/nf-core/busco/main'
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTQ_TRIM_FASTP_FASTQC     } from '../subworkflows/nf-core/fastq_trim_fastp_fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { SORTMERNA                   } from '../modules/nf-core/sortmerna/main'
include { SALMON_INDEX                } from '../modules/nf-core/salmon/index/main'
include { SALMON_QUANT                } from '../modules/nf-core/salmon/quant/main'
include { TRINITY                     } from '../modules/nf-core/trinity/main'
include { TRINITY as TRINITY_NO_NORM  } from '../modules/nf-core/trinity/main'


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
    // MODULE: FASTQ_TRIM_FASTP_FASTQC
    //
    FASTQ_TRIM_FASTP_FASTQC (
        INPUT_CHECK.out.reads,
        params.adapter_fasta,
        params.save_trimmed_fail,
        params.save_merged,
        params.skip_fastp,
        params.skip_fastqc
    )
    ch_filtered_reads = FASTQ_TRIM_FASTP_FASTQC.out.reads
    ch_versions = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions)

    if (!params.QC_only) {

        if (!params.quant_only) {
            //
            // MODULE: SORTMERNA
            //
            if (params.remove_ribo_rna) {
                ch_sortmerna_fastas = Channel.from(ch_ribo_db.readLines()).map { row -> file(row, checkIfExists: true) }.collect()

                SORTMERNA (
                    ch_filtered_reads,
                    ch_sortmerna_fastas
                )
                .reads
                .set { ch_filtered_reads }

                ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())
            }

            // All methods used pooled reads
            pool_ch = ch_filtered_reads.collect { meta, fastq -> fastq }.map { [[id:'pooled_reads', single_end:false], it] }      

            //
            // MODULE: CAT_FASTQ
            //
            CAT_FASTQ (
                pool_ch
            )
            ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

            // Method 1, 2 and 3 only use Trinity for assembly
            // Method 4 uses Trinity and rnaSPAdes for assembly
            // Method 5 combines assemblies from all methods
            // First Trinity assembly is created using normalized reads
            //
            // MODULE: Trinity 
            //
            TRINITY (
                CAT_FASTQ.out.reads
            )
            ch_versions = ch_versions.mix(TRINITY.out.versions)

            if (params.method == 3 | params.method == 5) {
                //
                // MODULE: Trinity (--no_normalize_reads)
                //
                TRINITY_NO_NORM (
                    CAT_FASTQ.out.reads
                )
            }  
            
            if (params.method == 4 | params.method == 5) {
                //
                // MODULE: RNASPADES
                //
                RNASPADES (
                    CAT_FASTQ.out.reads
                )
                ch_versions = ch_versions.mix(RNASPADES.out.versions)
            } 

            if (params.method == 1) {
                // Method 1 uses the trinity assembly as the final assembly 

                final_assembly_ch = TRINITY.out.trinity_assembly
                final_assembly_file = TRINITY.out.trinity_assembly.map{ meta, fasta -> fasta }

            } else if (params.method == 2) {
                // Method 2 does redundancy reduction on the trinity assembly
                // to create the final assembly
                //
                // MODULE: Evidential Gene
                //
                TR2AACDS (
                    TRINITY.out.trinity_assembly
                )
                ch_versions = ch_versions.mix(TR2AACDS.out.versions)

                final_assembly_ch = TR2AACDS.out.non_redundant_fasta
                final_assembly_file = TR2AACDS.out.non_redundant_fasta.map{ meta, fasta -> fasta }

            } else if (params.method == 3) {
                // Method 3 combines the trinity assembly with the trinity (--no_normalize) assembly
                // and does redundancy reduction on the combined assembly

                method_3_assemblies = TRINITY_NO_NORM.out.trinity_assembly
                    .mix(TRINITY.out.trinity_assembly)
                    .collect { meta, fasta -> fasta }
                    .map {[ [id:'all_assembled', single_end:false], it ] }

                //
                // MODULE: Evidential Gene
                //
                TR2AACDS (
                    method_3_assemblies
                )
                ch_versions = ch_versions.mix(TR2AACDS.out.versions)
                
                final_assembly_ch = TR2AACDS.out.non_redundant_fasta
                final_assembly_file = TR2AACDS.out.non_redundant_fasta.map{ meta, fasta -> fasta }

            } else if (params.method == 4) {
                // Method 4 combines the trinity assembly with the rnaSPAdes assembly (transcripts)
                // and does redundancy reduction on the combined assembly

                method_4_assemblies = RNASPADES.out.transcripts
                    .mix(TRINITY.out.trinity_assembly)
                    .collect { meta, fasta -> fasta }
                    .map {[ [id:'all_assembled', single_end:false], it ] }

                //
                // MODULE: Evidential Gene
                //
                TR2AACDS (
                    method_4_assemblies
                )
                ch_versions = ch_versions.mix(TR2AACDS.out.versions)
                
                final_assembly_ch = TR2AACDS.out.non_redundant_fasta
                final_assembly_file = TR2AACDS.out.non_redundant_fasta.map{ meta, fasta -> fasta }
                    
            } else if (params.method == 5) {
                // Method 5 combines the trinity assembly, the trinity (--no_normalize) assembly and the rnaSPAdes assembly
                // and does redundancy reduction on the combined assembly

                method_5_assemblies = TRINITY.out.trinity_assembly
                .mix(TRINITY_NO_NORM.out.trinity_assembly)
                .mix(RNASPADES.out.transcripts)
                .mix(RNASPADES.out.soft_filtered_transcripts)
                .mix(RNASPADES.out.hard_filtered_transcripts)
                .collect { meta, fasta -> fasta }
                .map {[ [id:'all_assembled', single_end:false], it ] }

                //
                // MODULE: Evidential Gene
                //
                TR2AACDS (
                    method_5_assemblies
                )
                ch_versions = ch_versions.mix(TR2AACDS.out.versions)
                
                final_assembly_ch = TR2AACDS.out.non_redundant_fasta
                final_assembly_file = TR2AACDS.out.non_redundant_fasta.map{ meta, fasta -> fasta }
            }

            //
            // MODULE: BUSCO
            //
            BUSCO (
                final_assembly_ch,
                "transcriptome",
                params.busco_lineage,
                params.busco_lineages_path,
                params.busco_config
            )
            ch_versions = ch_versions.mix(BUSCO.out.versions)

            //
            // MODULE: RNAQUAST
            //
            RNAQUAST (
                final_assembly_ch
            )
            ch_versions = ch_versions.mix(RNAQUAST.out.versions)
        }

        if (params.quant_only) {
            final_assembly_file = file(params.transcript_fasta, checkIfExists: true)
        }
        //
        // MODULE: Salmon index
        //
        SALMON_INDEX (
            final_assembly_file
        )
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)

        //
        // MODULE: Salmon quant
        //
        SALMON_QUANT (
            ch_filtered_reads,
            SALMON_INDEX.out.index,
            final_assembly_file,
            params.lib_type        
        )
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)

        //
        // MODULE: MERGE_TPM
        //
        MERGE_TPM (
            SALMON_QUANT.out.results.collect{it[1]}
        )
        ch_versions = ch_versions.mix(MERGE_TPM.out.versions)

        CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )
    }

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
    
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TRIM_FASTP_FASTQC.out.fastqc_raw_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TRIM_FASTP_FASTQC.out.fastqc_trim_zip.collect{it[1]}.ifEmpty([]))

    if (!params.QC_only) {
        if (params.remove_ribo_rna) {
            ch_multiqc_files = ch_multiqc_files.mix(SORTMERNA.out.log.collect{it[1]}.ifEmpty([]))
        }
        if (!params.quant_only) {
            ch_multiqc_files = ch_multiqc_files.mix(BUSCO.out.short_summaries_txt.collect{it[1]}.ifEmpty([]))
        }
        ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.results.collect{it[1]}.ifEmpty([]))
    }

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
