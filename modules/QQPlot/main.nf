process QQPlot {
    publishDir "${params.outdir}/plots/", mode: 'copy'

    input:
        path(assoc_tsv)
        val(pcol)

    output:
        path(qqplot), optional: true

    script:
        qqplot = "QQplot_${assoc_tsv.baseName}_${pcol}.png"

        """
        set -eo pipefail

        QQplot.R -i ${assoc_tsv} -p ${pcol} -o ${qqplot}
        """
}