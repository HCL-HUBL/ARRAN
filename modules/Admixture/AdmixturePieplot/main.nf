process AdmixturePieplot {
    publishDir "${params.outdir}/plots/", mode: 'copy'

    input:
        path(admixture_table)
        path(eigenvec)

    output:
        path(admixture_pieplot)

    script:
        admixture_pieplot = "${admixture_table}_pieplot.pdf"

        """
        set -eo pipefail

        admixturePie.R \
            -e ${eigenvec} \
            -a ${admixture_table} \
            -o ${admixture_pieplot}
        """
}