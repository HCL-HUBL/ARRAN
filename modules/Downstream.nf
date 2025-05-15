#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process ManhattanPlot {
    publishDir "${params.outdir}/plots/", mode: 'copy'

    input:
        path(saige_sv)

    output:
        path(manhattan)
    
    script:
        manhattan = "Manhattan_${saige_sv.baseName}.pdf"
        
        """
        set -eo pipefail

        ${params.tools.Rscript} ${projectDir}/bin/manhattan_plot.R -i ${saige_sv} -o ${manhattan}
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
