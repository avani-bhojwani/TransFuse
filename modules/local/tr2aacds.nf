process TR2AACDS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::perl=5.34.0 bioconda::exonerate=2.4.0 bioconda::cd-hit=4.8.1 bioconda::blast=2.12.0"
    container "$launchDir/containers/evigene.sif"

    input:
    tuple val(meta), path(assemblies)

    output:
    tuple val(meta), path("*.okay.fa")  , emit: non_redundant_fasta
    path "versions.yml"                 , emit: versions

    script:
    def mem_MB = (task.memory.toMega())
    """
    cat ${assemblies} > ${meta.id}.fa
    gzip ${meta.id}.fa

    tr2aacds.pl -NCPU $task.cpus -MAXMEM ${mem_MB} -log -cdna ${meta.id}.fa.gz

    cp okayset/*.okay.mrna ${meta.id}.okay.fa

    # replace headers with the old ones
    sed -n -i "/^>/s/.*oid=\\([^;]*\\);/>\\1/p; t; p" ${meta.id}.okay.fa

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