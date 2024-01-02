# ![nf-core/transfuse](docs/images/nf-core-transfuse_logo_light.png#gh-light-mode-only) ![nf-core/transfuse](docs/images/nf-core-transfuse_logo_dark.png#gh-dark-mode-only)

[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/transfuse/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/transfuse)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23transfuse-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/transfuse)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/transfuse** is a bioinformatics pipeline for de novo transcriptome assembly of paired-end short reads. It takes a samplesheet and FASTQ files as input, perfoms quality control (QC), trimming, assembly, redundancy reduction, pseudoalignment, and quantification. It outputs a transcriptome assembly FASTA file, a transcript abundance TSV file, and a MultiQC report with assembly quality and read QC metrics.

![nf-core/transfuse metro map](docs/images/transfuse_metro_map.drawio.svg)

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Read QC of raw reads ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter and quality trimming ([`fastp`](https://github.com/OpenGene/fastp))
3. Remove rRNA or mitochondrial DNA (optional) ([`SortMeRNA`](https://hpc.nih.gov/apps/sortmeRNA.html))
4. Read QC of trimmed reads ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
5. Transcriptome assembly and redundancy reduction (choose a method):
   1. Assembly with [`Trinity`](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
   2. Assembly with [`Trinity`](https://github.com/trinityrnaseq/trinityrnaseq/wiki), redundancy reduction with [`Evidential Gene tr2aacds`](http://arthropods.eugenes.org/EvidentialGene/) (default)
   3. Assembly with [`Trinity`](https://github.com/trinityrnaseq/trinityrnaseq/wiki) and [`rnaSPAdes`](https://cab.spbu.ru/software/rnaspades/) (with multiple k-mer sizes), redundancy reduction with [`Evidential Gene tr2aacds`](http://arthropods.eugenes.org/EvidentialGene/)
6. Assembly completeness QC ([`BUSCO`](https://busco.ezlab.org/))
7. Other assembly quality metrics ([`rnaQUAST`](http://cab.spbu.ru/software/rnaquast))
8. Pseudo-alignment and quantification ([`Salmon`](https://combine-lab.github.io/salmon/))
9. Present HTML report for raw reads, trimmed reads, BUSCO, and Salmon ([`MultiQC`](http://multiqc.info/))

## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):
-->
First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,SRX9161249_SRR12681131_1.fastq.gz,SRX9161255_SRR12681131_2.fastq.gz
CONTROL_REP2,SRX9161250_SRR12681132_1.fastq.gz,SRX9161251_SRR12681132_2.fastq.gz
CONTROL_REP3,SRX9161251_SRR12681133_1.fastq.gz,SRX9161251_SRR12681133_2.fastq.gz
```

Each row represents a pair of fastq files (paired end).

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run nf-core/transfuse \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details, please refer to the [usage documentation](https://nf-co.re/transfuse/usage) and the [parameter documentation](https://nf-co.re/transfuse/parameters).

## Pipeline output

To see the the results of a test run with a full size dataset refer to the [results](https://nf-co.re/transfuse/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/transfuse/output).

## Credits

nf-core/transfuse was originally written by Avani Bhojwani.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#transfuse` channel](https://nfcore.slack.com/channels/transfuse) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/transfuse for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
