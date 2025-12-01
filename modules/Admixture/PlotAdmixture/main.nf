process PlotAdmixture {
    publishDir "${params.outdir}/plots/", mode: 'copy'

    input:
        path(admixture_table)

    output:
        path(admixture_barplot)

    script:
        admixture_barplot = "${admixture_table}.pdf"

        """
        set -eo pipefail

        admixture.R \
            -i ${admixture_table} \
            -k ${params.admixture_K} \
            -o ${admixture_barplot}
        """
}