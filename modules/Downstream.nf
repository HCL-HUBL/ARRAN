#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process ManhattanPlot {
    publishDir "${params.outdir}/plots/", mode: 'copy'

    input:
        path(assoc_tsv)
        val(chr_col)
        val(pos_col)
        val(marker_col)
        val(pval_col)

    output:
        path(manhattan)
    
    script:
        manhattan = "Manhattan_${assoc_tsv.baseName}.pdf"
        
        """
        set -eo pipefail

        ${params.tools.Rscript} ${projectDir}/bin/manhattan_plot.R -i ${assoc_tsv} \
            -c ${chr_col} -b ${pos_col} -m ${marker_col} -p ${pval_col} \
            -o ${manhattan}
        """
}


process QQPlot {
    publishDir "${params.outdir}/plots/", mode: 'copy'

    input:
        path(assoc_tsv)
        val(pcol)

    output:
        path(qqplot)

    script:
        qqplot = "QQplot_${assoc_tsv.baseName}.png"

        """
        set -eo pipefail

        ${params.tools.Rscript} ${projectDir}/bin/QQplot.R -i ${assoc_tsv} -p ${pcol} -o ${qqplot}
        """
}

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

        ${params.tools.Rscript} ${projectDir}/bin/create_summary_stats.R -i ${assoc_tsv} -s ${software} -t ${params.trait_type} -o ${summ_stats}
        """
}