process AdmixtureBarplot {
    publishDir "${params.outdir}/plots/", mode: 'copy'

    label 'admixture'
    cpus 1
    memory { 100.MB * task.attempt }
    time { 1.minute * task.attempt }

    input:
        path(admixture_table)

    output:
        path(admixture_barplot)

    script:
        admixture_barplot = "${admixture_table}_barplot.pdf"

        """
        set -eo pipefail

        admixture.R \
            -i ${admixture_table} \
            -k ${params.admixture_K} \
            -o ${admixture_barplot}
        """
}