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
include { RNAQUAST                    } from '../modules/local/rnaquast'
include { ASSEMBLE; ASSEMBLE as ASSEMBLE_FIRST_SAMPLE } from '../subworkflows/local/assemble'

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
include { BUSCO                       } from '../modules/nf-core/busco/main'
include { TR2AACDS; TR2AACDS as FIRST_TR2AACDS } from '../modules/local/tr2aacds'
include { CAT_FASTA                   } from '../modules/local/cat_fasta'

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

    // Method 1, method 4, and method 5 pool all reads together and then assemble them
    if (params.method == 1 | params.method == 4 | params.method == 5) {
        method_1_pool_ch = FASTP.out.reads.collect { meta, fastq -> fastq }.map { [[id:'pooled_reads', single_end:false], it] }      

        //
        // MODULE: CAT_FASTQ
        //
        CAT_FASTQ (
            method_1_pool_ch
        )
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

        // Method 4 and method 5 only use Trinity for assembly
        if (params.method == 4 | params.method == 5) {
            //
            // MODULE: Trinity
            //
            TRINITY (
                CAT_FASTQ.out.reads
            )
            ch_versions = ch_versions.mix(TRINITY.out.versions)

            // Method 4 uses the trinity assembly as the final assembly
            if (params.method == 4) {
                final_assembly_ch = TRINITY.out.trinity_assembly
                final_assembly_file = TRINITY.out.trinity_assembly.map{ meta, fasta -> fasta }
            } else if (params.method == 5) {
                // Method 5 does redundancy reduction on the trinity assembly
                //
                // MODULE: Evidential GeneF
                //
                TR2AACDS (
                    TRINITY.out.trinity_assembly
                )
                ch_versions = ch_versions.mix(TR2AACDS.out.versions)

                final_assembly_ch = TR2AACDS.out.non_redundant_fasta
                final_assembly_file = TR2AACDS.out.non_redundant_fasta.map{ meta, fasta -> fasta }
            }
        
        // Method 1 uses Trinity, rnaSPAdes, and evigene's tr2aacds for assembly and redundancy reduction
        } else if (params.method == 1) {
            //
            // MODULE: ASSEMBLE
            //
            ASSEMBLE (
                CAT_FASTQ.out.reads,
                params.kmers
            )
            ch_versions = ch_versions.mix(ASSEMBLE.out.versions)

            method_1_assemblies = ASSEMBLE.out.trinity_assembly.mix(ASSEMBLE.out.spades_assembly).collect { meta, fasta -> fasta }.map {[ [id:'all_assembled', single_end:false], it ] }

            //
            // MODULE: Evidential Gene
            //
            TR2AACDS (
                method_1_assemblies
            )
            ch_versions = ch_versions.mix(TR2AACDS.out.versions)
            
            final_assembly_ch = TR2AACDS.out.non_redundant_fasta
            final_assembly_file = TR2AACDS.out.non_redundant_fasta.map{ meta, fasta -> fasta }
        }

    } else if (params.method == 2) {
        // Method 2 assembles each sample separately, them combines the assemblies
        //
        // MODULE: ASSEMBLE
        //
        ASSEMBLE (
            FASTP.out.reads,
            params.kmers
        )
        ch_versions = ch_versions.mix(ASSEMBLE.out.versions)

        method_2_assemblies = ASSEMBLE.out.trinity_assembly.mix(ASSEMBLE.out.spades_assembly)
        method_2_assemblies = method_2_assemblies.collect { meta, fasta -> fasta }.map {[ [id:'all_assembled', single_end:false], it ] }

        //
        // MODULE: Evidential Gene
        //
        TR2AACDS (
            method_2_assemblies
        )
        ch_versions = ch_versions.mix(TR2AACDS.out.versions)

        final_assembly_ch = TR2AACDS.out.non_redundant_fasta
        final_assembly_file = TR2AACDS.out.non_redundant_fasta.map{ meta, fasta -> fasta }

    } else if ( params.method == 3 ) {
        // Method 3 creates a 'reference' transcriptome from one sample
        first_sample_ch = FASTP.out.reads.filter{ it[0].first == true }
        remaining_samples_ch = FASTP.out.reads.filter{ it[0].first == false }

        //
        // MODULE: ASSEMBLE FIRST SAMPLE
        //
        ASSEMBLE_FIRST_SAMPLE (
            first_sample_ch,
            params.kmers
        )
        ch_versions = ch_versions.mix(ASSEMBLE_FIRST_SAMPLE.out.versions)

        first_sample_assemblies = ASSEMBLE_FIRST_SAMPLE.out.trinity_assembly.mix(ASSEMBLE_FIRST_SAMPLE.out.spades_assembly)
        first_sample_assemblies = first_sample_assemblies.collect { meta, fasta -> fasta }.map {[ [id:'first_sample', single_end:false], it ] }
        //
        // MODULE: Evidential Gene for first sample assembly
        //
        FIRST_TR2AACDS (
            first_sample_assemblies
        )
        ch_versions = ch_versions.mix(FIRST_TR2AACDS.out.versions)

        first_sample_assembly_ch = FIRST_TR2AACDS.out.non_redundant_fasta
        first_sample_assembly_file = FIRST_TR2AACDS.out.non_redundant_fasta.map{ meta, fasta -> fasta }

        // Align all remaining reads to the reference transcriptome
        // 
        // MODULE: STAR genomeGenerate
        //
        STAR_GENOMEGENERATE (
            first_sample_assembly_file,
            params.gtf,
            params.star_ignore_sjdbgtf
        )
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

        //
        // MODULE: STAR align
        //
        STAR_ALIGN (
            remaining_samples_ch,
            STAR_GENOMEGENERATE.out.index,
            params.gtf,
            params.star_ignore_sjdbgtf,
            params.seq_platform,
            params.seq_center
        )
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

        // Pool unmapped reads from remaining samples before assembly
        method_3_pool_ch = STAR_ALIGN.out.fastq.collect { meta, fastq -> fastq }.map {[ [id:'pooled_unmapped_reads', single_end:false], it ] }

        //
        // MODULE: CAT_FASTQ for remaining samples
        //
        CAT_FASTQ (
            method_3_pool_ch
        )

        //
        // MODULE: ASSEMBLE remaining samples
        //
        ASSEMBLE (
            CAT_FASTQ.out.reads,
            params.kmers
        )
        ch_versions = ch_versions.mix(ASSEMBLE.out.versions)

        method_3_remaining_assemblies = ASSEMBLE.out.trinity_assembly.mix(ASSEMBLE.out.spades_assembly).collect { meta, fasta -> fasta }.map {[ [id:'remaining_samples', single_end:false], it ] }
        
        //
        // MODULE: Evidential Gene for remaining samples
        //
        TR2AACDS (
            method_3_remaining_assemblies
        )
        ch_versions = ch_versions.mix(TR2AACDS.out.versions)

        all_assemblies_ch = TR2AACDS.out.non_redundant_fasta.mix(first_sample_assembly_ch).collect { meta, fasta -> fasta }.map {[ [id:'all_assembled', single_end:false], it ] }

        //
        // Concat all assemblies together
        //
        CAT_FASTA (
            all_assemblies_ch
        )

        final_assembly_ch = CAT_FASTA.out.merged_assembly
        final_assembly_file = CAT_FASTA.out.merged_assembly.map{ meta, fasta -> fasta }
    } 
    //
    // MODULE: BUSCO
    //
    BUSCO (
        final_assembly_ch,
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
        FASTP.out.reads,
        SALMON_INDEX.out.index,
        final_assembly_file,
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

    if ( params.method == 3 ) {
        ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect{it[1]}.ifEmpty([]))
    }
    ch_multiqc_files = ch_multiqc_files.mix(BUSCO.out.short_summaries_txt.collect{it[1]}.ifEmpty([]))

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
