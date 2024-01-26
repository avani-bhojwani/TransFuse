# nf-core/transfuse Pipeline Parameters Documentation

A pipeline for de novo transcriptome assembly.

## Input/output options

Define where the pipeline should find input data and save output data.

- **input**:
    - **Type**: string
    - **Format**: file-path
    - **Mimetype**: text/csv
    - **Pattern**: `^\S+\.csv$`
    - **Description**: Path to comma-separated file containing information about the samples in the experiment. See README for example.

- **outdir**:
    - **Type**: string
    - **Format**: directory-path
    - **Description**: The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.

- **email**:
    - **Type**: string
    - **Description**: Email address for completion summary.

- **multiqc_title**:
    - **Type**: string
    - **Description**: MultiQC report title. Printed as page header, used for filename if not otherwise specified.


## FastQC/FastP options

- **adapter_fasta**:
    - **Type**: string
    - **Default**: `[]`
    - **Format**: file-path
    - **Description**: File in FASTA format containing possible adapters to remove. Accepted formats: *.{fasta,fna,fas,fa}.

- **save_trimmed_fail**:
    - **Type**: boolean
    - **Description**: Specify true to save files that failed to pass trimming thresholds ending in *.fail.fastq.gz.

- **save_merged**:
    - **Type**: boolean
    - **Description**: Specify true to save all merged reads to a file ending in *.merged.fastq.gz.

- **skip_fastp**:
    - **Type**: boolean
    - **Description**: Skip the fastp process if true.

- **skip_fastqc**:
    - **Type**: boolean
    - **Description**: Skip FastQC processes if true.

- **extra_fastp_args**:
    - **Type**: string
    - **Description**: Extra arguments for fastp. For example, '--trim_front1 15 --trim_front2 15 --trim_tail1 5 --trim_tail2 5'.

- **QC_only**:
    - **Type**: boolean
    - **Description**: Only run FastQC and fastp steps. Assembly and subsequent steps won't be performed.


## SortMeRNA options

Options related to SortMeRNA for ribosomal RNA removal.

- **remove_ribo_rna**:
    - **Type**: boolean
    - **Description**: Enable the removal of reads derived from ribosomal RNA using SortMeRNA.

- **save_non_ribo_reads**:
    - **Type**: boolean
    - **Description**: If this option is specified, intermediate FastQ files containing non-rRNA reads will be saved in the results directory.

- **ribo_database_manifest**:
    - **Type**: string
    - **Default**: `${projectDir}/assets/rrna-db-defaults.txt`
    - **Description**: Text file containing paths to fasta files (one per line) that will be used to create the database for SortMeRNA.

## Assembly options

Options related to the assembly process.

- **method**:
    - **Type**: integer
    - **Default**: 2
    - **Description**: Method to use for assembly. See README for details.

- **extra_trinity_args**:
    - **Type**: string
    - **Description**: Extra arguments to pass to Trinity command in addition to defaults.

- **extra_tr2aacds_args**:
    - **Type**: string
    - **Description**: Extra arguments for tr2aacds.pl. For example, '-MINAA=20'.

## BUSCO options

Options related to BUSCO (Benchmarking Universal Single-Copy Orthologs) analysis.

- **busco_lineage**:
    - **Type**: string
    - **Default**: `arthropoda_odb10.2019-11-20`
    - **Description**: The BUSCO lineage to use for the analysis, or `"auto"` to automatically select the lineage.

- **busco_lineages_path**:
    - **Type**: string
    - **Default**: `[]`
    - **Description**: Path to local BUSCO lineages directory.
    - **Format**: directory-path

- **busco_config**:
    - **Type**: string
    - **Default**: `[]`
    - **Description**: Path to BUSCO config file.
    - **Format**: file-path

## Salmon options

Options related to the Salmon tool for transcript quantification.

- **lib_type**:
    - **Type**: string
    - **Default**: `A`
    - **Description**: Override library type inferred based on strandedness defined in the meta object. 

## Institutional config options

Parameters used to describe centralized config profiles. These should not be edited.

**Help Text**: The centralized nf-core configuration profiles use a handful of pipeline parameters to describe themselves. This information is then printed to the Nextflow log when you run a pipeline. You should not need to change these values when you run a pipeline.

- **custom_config_version**:
    - **Type**: string
    - **Default**: `master`
    - **Description**: Git commit id for Institutional configs.
    - **Hidden**: true

- **custom_config_base**:
    - **Type**: string
    - **Default**: `https://raw.githubusercontent.com/nf-core/configs/master`
    - **Description**: Base directory for Institutional configs.
    - **Hidden**: true
    - **Help Text**: If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.

- **config_profile_name**:
    - **Type**: string
    - **Description**: Institutional config name.
    - **Hidden**: true

- **config_profile_description**:
    - **Type**: string
    - **Description**: Institutional config description.
    - **Hidden**: true

- **config_profile_contact**:
    - **Type**: string
    - **Description**: Institutional config contact information.
    - **Hidden**: true

- **config_profile_url**:
    - **Type**: string
    - **Description**: Institutional config URL link.
    - **Hidden**: true

## Max job request options

Set the top limit for requested resources for any single job.

**Help Text**: If you are running on a smaller system, a pipeline step requesting more resources than are available may cause Nextflow to stop the run with an error. These options allow you to cap the maximum resources requested by any single job so that the pipeline will run on your system.

Note that you can not _increase_ the resources requested by any job using these options. For that, you will need your own configuration file. See [the nf-core website](https://nf-co.re/usage/configuration) for details.

- **max_cpus**:
    - **Type**: integer
    - **Default**: 16
    - **Description**: Maximum number of CPUs that can be requested for any single job.
    - **Hidden**: true
    - **Help Text**: Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`.

- **max_memory**:
    - **Type**: string
    - **Default**: `128.GB`
    - **Description**: Maximum amount of memory that can be requested for any single job.
    - **Pattern**: `^\d+(\.\d+)?\.?\s*(K|M|G|T)?B$`
    - **Hidden**: true
    - **Help Text**: Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`.

- **max_time**:
    - **Type**: string
    - **Default**: `240.h`
    - **Description**: Maximum amount of time that can be requested for any single job.
    - **Pattern**: `^(\d+\.?\s*(s|m|h|day)\s*)+$`
    - **Hidden**: true
    - **Help Text**: Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`.

## Generic options

Less common options for the pipeline, typically set in a config file.

**Help Text**: These options are common to all nf-core pipelines and allow you to customize some of the core preferences for how the pipeline runs.

Typically, these options would be set in a Nextflow config file loaded for all pipeline runs, such as `~/.nextflow/config`.

- **help**:
    - **Type**: boolean
    - **Description**: Display help text.

- **version**:
    - **Type**: boolean
    - **Description**: Display version and exit.

- **publish_dir_mode**:
    - **Type**: string
    - **Default**: `copy`
    - **Description**: Method used to save pipeline results to output directory.
    - **Help Text**: The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.
    - **Enum**: ["symlink", "rellink", "link", "copy", "copyNoFollow", "move"]

- **email_on_fail**:
    - **Type**: string
    - **Description**: Email address for completion summary, only when the pipeline fails.
    - **Help Text**: An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.
    - **Pattern**: `^([a-zA-Z0-9_\-\.]+)@([a-zA-Z0-9_\-\.]+)\.([a-zA-Z]{2,5})$`

- **plaintext_email**:
    - **Type**: boolean
    - **Description**: Send plain-text email instead of HTML.

- **max_multiqc_email_size**:
    - **Type**: string
    - **Default**: `25.MB`
    - **Description**: File size limit when attaching MultiQC reports to summary emails.
    - **Pattern**: `^\d+(\.\d+)?\.?\s*(K|M|G|T)?B$`

- **monochrome_logs**:
    - **Type**: boolean
    - **Description**: Do not use coloured log outputs.

- **hook_url**:
    - **Type**: string
    - **Description**: Incoming hook URL for messaging service.
    - **Help Text**: Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.

- **multiqc_config**:
    - **Type**: string
    - **Description**: Custom config file to supply to MultiQC.

- **multiqc_logo**:
    - **Type**: string
    - **Description**: Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file.

- **multiqc_methods_description**:
    - **Type**: string
    - **Description**: Custom MultiQC yaml file containing HTML including a methods description.

- **tracedir**:
    - **Type**: string
    - **Default**: `${params.outdir}/pipeline_info`
    - **Description**: Directory to keep pipeline Nextflow logs and reports.

- **validate_params**:
    - **Type**: boolean
    - **Default**: true
    - **Description**: Boolean whether to validate parameters against the schema at runtime.

- **show_hidden_params**:
    - **Type**: boolean
    - **Description**: Show all params when using `--help`.
    - **Help Text**: By default, parameters set as _hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters.
