process QQPlot {
    publishDir "${params.outdir}/plots/", mode: 'copy'

    label 'GWAS'
    cpus 1
    memory { 1.GB * task.attempt }
    time { 5.minute * task.attempt }
    
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