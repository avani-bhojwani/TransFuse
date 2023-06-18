process RNASPADES {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::spades=3.15.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spades:3.15.5--h95f258a_1' :
        'biocontainers/spades:3.15.5--h95f258a_1' }"

    input:
    val kmers
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.SPADES.fa")       , emit: spades_assembly
    path "versions.yml"                    , emit: versions

    script:
    """
    mem=\$( echo ${task.memory} | cut -f 1 -d " " )

    #run rnaSPAdes for each kmer, and stores the output in a folder 
    #the folder's name  includes the sample ID and the current k value
    for kmer in `echo $kmers | tr "," " "`;do
        rnaspades.py -1 ${reads[0]} -2 ${reads[1]} -o ${meta.id}_spades_\${kmer} -t ${task.cpus} -k \${kmer} -m \${mem}
    done

    #modify the identifiers of the sequences in the transcripts.fasta files, 
    #adding a prefix that specifies the k-mer size used in the assembly
    for kmer in `echo $kmers | tr "," " "`;do
        sed -i "s/>/>SPADES.k\${kmer}./g" ${meta.id}_spades_\${kmer}/transcripts.fasta
    done

    #concat all the fasta files from the individual assemblies into a single file
    cat ${meta.id}_spades_*/transcripts.fasta >${meta.id}.SPADES.fa

    #copy each fasta file to a new file that specifies the k-mer size in its name.
    for kmer in `echo $kmers | tr "," " "`;do
        cp ${meta.id}_spades_\${kmer}/transcripts.fasta ${meta.id}.SPADES.k\${kmer}.fa
    done

    #remove the individual assembly folders
    rm -r ${meta.id}_spades_*
    
    v=\$( rnaspades.py -v 2>&1 | awk '{print \$4}' | tr -d "v" )
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaSPAdes: \$v
    END_VERSIONS
    """
}