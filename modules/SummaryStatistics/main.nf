process SummaryStatistics {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "$it" }, mode: 'copy'

    label 'GWAS'
    cpus 1
    memory { 500.MB * task.attempt }
    time { 5.minute * task.attempt }
    
    input:
        path(assoc_tsv)
        val(software)       //Software that was used to generate the assoc_tsv file

    output:
        path(summ_stats)
    
    script:
        summ_stats = "Summary_stats_${assoc_tsv}"

        """
        set -eo pipefail

        create_summary_stats.R -i ${assoc_tsv} -s ${software} -t ${params.trait_type} -o ${summ_stats}
        """
}
