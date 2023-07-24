//
// Assemble clean reads using Trinity and rnaSPAdes
//

include { TRINITY                     } from '../../modules/nf-core/trinity/main'
include { RNASPADES                   } from '../../modules/local/rnaspades'

workflow ASSEMBLE {
    take:
    reads           // channel: [ val(meta), [ reads ] ]
    kmers           // string of numbers seprated by commas

    main:
    ch_versions    = Channel.empty()
    
    TRINITY ( reads ).trinity_assembly.set { trinity_assembly }
    ch_versions = ch_versions.mix(TRINITY.out.versions.first())

    RNASPADES ( kmers, reads ).spades_assembly.set { spades_assembly }
    ch_versions = ch_versions.mix(RNASPADES.out.versions)

    emit:
    trinity_assembly        // channel: [ val(meta), [ assembly ] ]
    spades_assembly         // channel: [ val(meta), [ assembly ] ]
    versions = ch_versions  // channel: [ versions.yml ]
}