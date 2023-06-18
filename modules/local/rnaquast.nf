process RNAQUAST {
    label 'process_medium'

    conda "bioconda::rnaquast=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rnaquast:2.2.3--h9ee0642_0' :
        'biocontainers/rnaquast:2.2.3--h9ee0642_0' }"

    input:
        path fasta

    output:
        path "logs"            , emit: logs
        path "*_output"         , emit: output
        path "short_report.pdf" , emit: report_pdf
        path "short_report.tex" , emit: report_tex
        path "short_report.tsv" , emit: report_tsv
        path "short_report.txt" , emit: report_txt
        path "versions.yml"     , emit: versions

    script:
    """
    rnaQUAST.py --transcripts $fasta -o .
    v=\$(rnaQUAST.py | grep "QUALITY ASSESSMENT" | head -n1 | awk -F " v." '{print \$2}')
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaQUAST: \$v
    END_VERSIONS
    """
}