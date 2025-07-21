#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// chrX-specific QC (only on controls, so only for binary traits):
// XWAS will filter variants which have significantly =/= MAFs between males and females in controls:
// 
process ChrX_specific_QC {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "xwas/$it" }, mode: 'copy'
    
    input:
        val(plink_chrX_basename)
        path(plink_chrX_bed)
        path(plink_chrX_bim)
        path(plink_chrX_fam)
        val(n_vars)

    output:
        tuple val(plink_chrX_QCed_basename), path(plink_chrX_QCed_files), emit: plink_chrX_QCed
        path("ChrX_specific_QC.log")

    script:
        plink_chrX_QCed_basename = "${params.out_basename}_chrX_QCed"
        plink_chrX_QCed_files    = "${params.out_basename}_chrX_QCed.{bim,bed,fam}"

        bonfx = params.xwas_alpha / n_vars

        """
        set -eo pipefail

        ${params.tools.xwas} --noweb --xwas \
            --bfile ${plink_chrX_basename} \
            --make-bed --out ${plink_chrX_QCed_basename} \
            --freqdiff-x ${bonfx} > ChrX_specific_QC.log
        """
}

process ChrX_SNVs_Assoc {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "xwas/$it" }, mode: 'copy'

    input:
        tuple val(plink_chrX_QCed_basename), path(plink_chrX_QCed_files)
        path(phenoFile)

    output:
        path(x_output), emit: xstrat_assoc

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

