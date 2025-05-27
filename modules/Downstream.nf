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
        path(saige_sv)
        val(pcol)

    output:
        path(qqplot)

    script:
        qqplot = "QQplot_${saige_sv.baseName}.png"

        """
        set -eo pipefail

        ${params.tools.Rscript} ${projectDir}/bin/QQplot.R -i ${saige_sv} -p ${pcol} -o ${qqplot}
        """
}
