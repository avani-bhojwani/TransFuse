process RNAQUAST {
    label 'process_medium'

    conda "bioconda::rnaquast=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rnaquast:2.2.3--h9ee0642_0' :
        'biocontainers/rnaquast:2.2.3--h9ee0642_0' }"

    input:
        path fasta

    output:
        path "rnaQUAST_results", emit: rnaquast_results
        path "versions.yml"    , emit: versions

    script:
    """
    rnaQUAST.py --transcripts $fasta
    v=\$(rnaQUAST.py | grep "QUALITY ASSESSMENT" | head -n1 | awk -F " v." '{print \$2}')
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaQUAST: \$v
    END_VERSIONS
    """
}