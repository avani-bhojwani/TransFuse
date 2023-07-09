process TRFORMAT {
    label 'process_single'

    conda "bioconda::perl=5.34.0 bioconda::exonerate=2.4.0 bioconda::cd-hit=4.8.1 bioconda::blast=2.12.0"
    container "$launchDir/containers/evigene.sif"

    input:
    path old_fasta
    tuple val(meta), path(trinity_assembly)
    tuple val(meta), path(spades_assembly)

    output:
    path "*.tr", emit: reformatted_fasta
    path "versions.yml", emit: versions

    script:
    """
    gunzip $trinity_assembly -c > trinity.fa

    # when mapping to reference transcriptome 
    if [ !${params.skip_mapping} ]; then  
        trformat.pl -pre $task.process -output all.tr -input $old_fasta trinity.fa $spades_assembly
    #when skipping mapping to reference transcriptome
    else
        trformat.pl -pre $task.process -output all.tr -input trinity.fa $spades_assembly
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl --version | head -2 | tail -1 | cut -d'(' -f2 | cut -d')' -f1)
        exonerate: \$(exonerate --version | grep exonerate | cut -d' ' -f5)
        cd-hit: \$(cd-hit | grep "CD-HIT version" | cut -d' ' -f4)
        blast: \$(blastn -version | grep -i 'blastn: ' | cut -d' ' -f2)
        EvidentialGene: 'accessed 26/05/2023'
    END_VERSIONS
    """
}