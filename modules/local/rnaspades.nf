process RNASPADES {
    tag "$meta.id"
    label 'process_high_memory'

    conda "bioconda::spades=3.15.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spades:3.15.5--h95f258a_1' :
        'biocontainers/spades:3.15.5--h95f258a_1' }"

    input:
    val kmers
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*transcripts.fasta")       , emit: spades_assembly
    path "versions.yml"                    , emit: versions

    script:
    """
    mem=\$( echo ${task.memory} | cut -f 1 -d " " )
    sample_id=${meta.id}

    rnaspades.py -1 ${reads[0]} -2 ${reads[1]} -o ${meta.id}_spades -t ${task.cpus} -m \${mem}

    cp ${meta.id}_spades/transcripts.fasta ${meta.id}_transcripts.fasta
    
    v=\$( rnaspades.py -v 2>&1 | awk '{print \$4}' | tr -d "v" )
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaSPAdes: \$v
    END_VERSIONS
    """
}