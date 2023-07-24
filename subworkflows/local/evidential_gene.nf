//
// Reduce redundancy using Evidential Gene
//

include { TRFORMAT                    } from '../../modules/local/trformat'
include { TR2AACDS                    } from '../../modules/local/tr2aacds'

workflow EVIDENTIAL_GENE {
    take:
    assemblies              // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions    = Channel.empty()

    TRFORMAT ( assemblies )
    TR2AACDS ( TRFORMAT.out.reformatted_fasta )
    ch_versions = ch_versions.mix(TR2AACDS.out.versions)

    emit:
    non_redundant_fasta = TR2AACDS.out.non_redundant_fasta
    versions = ch_versions  // channel: [ versions.yml ]
}