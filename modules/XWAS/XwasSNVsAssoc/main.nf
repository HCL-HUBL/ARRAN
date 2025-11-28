process XwasSNVsAssoc {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "xwas/$it" }, mode: 'copy'

    input:
        tuple val(plink_chrX_QCed_basename), path(plink_chrX_QCed_files)
        path(phenoFile)

    output:
        path(x_output), emit: xstrat_assoc
        path("xwas.log")

    script:
        if(params.trait_type == "binary") {
            x_output = "xwas.xstrat.logistic"
            cmd_beta = ""
        } else {
            x_output = "xwas.xstrat.linear"
            cmd_beta = "--xbeta" // Transform Odd Ratios into betas in case of a quantitative traits
        }

        """
        set -eo pipefail

        # --xchr-model 2 will code male genotypes as 0/2
        # using --fishers method to combine p-values:
        ${params.tools.xwas} --noweb --xwas \
            --bfile ${plink_chrX_QCed_basename} \
            --covar ${phenoFile} --covar-name ${params.xwas_covar} \
            --strat-sex --fishers --xchr-model 2 --ci 0.95 \
            ${cmd_beta} 
        """
}