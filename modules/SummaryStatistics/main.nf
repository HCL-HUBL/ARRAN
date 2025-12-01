process SummaryStatistics {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "$it" }, mode: 'copy'

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