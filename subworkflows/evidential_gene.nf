//
// Reduce redundancy using Evidential Gene
//

include { TRFORMAT                    } from '../modules/local/trformat'
include { TR2AACDS                    } from '../modules/local/tr2aacds'

workflow EVIDENTIAL_GENE {
    take:
    assemblies              // channel: [ val(meta), [ fasta ] ]

    main:
    TRFORMAT ( assemblies )
    TR2AACDS ( TRFORMAT.reformatted_fasta ).non_redundant_fasta.set { non_redundant_fasta }
    ch_versions = ch_versions.mix(TR2AACDS.out.versions)

    emit:
    non_redundant_fasta     // channel: [ non_redundant.okay.fa ]
    versions = ch_versions  // channel: [ versions.yml ]
}