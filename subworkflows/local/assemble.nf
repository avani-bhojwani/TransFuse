//
// Assemble clean reads using Trinity and rnaSPAdes
//

include { TRINITY                     } from '../modules/nf-core/trinity/main'
include { RNASPADES                   } from '../modules/local/rnaspades'
include { TRFORMAT                    } from '../modules/local/trformat'
include { TR2AACDS                    } from '../modules/local/tr2aacds'

workflow ASSEMBLE {
    take:
    reads           // channel: [ val(meta), [ reads ] ]
    kmers           // string of numbers seprated by commas

    main:
    TRINITY ( reads ).trinity_assembly.set { trinity_assembly }
    ch_versions = ch_versions.mix(TRINITY.out.versions.first())

    RNASPADES ( kmers, reads ).spades_assembly.set { spades_assembly }
    ch_version = ch_versions.mix(RNASPADES.out.versions)

    TRFORMAT ( trinity_assembly, spades_assembly )
    TR2AACDS ( TRFORMAT.reformatted_fasta ).non_redundant_fasta.set { non_redundant_fasta }
    ch_versions = ch_versions.mix(TR2AACDS.out.versions)

    emit:
    trinity_assembly        // channel: [ val(meta), [ assembly ] ]
    spades_assembly         // channel: [ val(meta), [ assembly ] ]
    non_redundant_fasta     // channel: [ non_redundant.okay.fa ]
    versions = ch_versions  // channel: [ versions.yml ]
}