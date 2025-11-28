// chrX-specific QC (only on controls, so only for binary traits):
// XWAS will filter variants which have significantly =/= MAFs between males and females in controls:
process XwasQC {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "xwas/$it" }, mode: 'copy'
    
    input:
        val(plink_chrX_basename)
        path(plink_chrX_bed)
        path(plink_chrX_bim)
        path(plink_chrX_fam)
        val(n_vars)

    output:
        tuple val(plink_chrX_QCed_basename), path(plink_chrX_QCed_files), emit: plink_chrX_QCed
        path("XwasQC.log")

    script:
        plink_chrX_QCed_basename = "${params.out_basename}_chrX_QCed"
        plink_chrX_QCed_files    = "${params.out_basename}_chrX_QCed.{bim,bed,fam}"

        bonfx = params.xwas_alpha / n_vars

        """
        set -eo pipefail

        ${params.tools.xwas} --noweb --xwas \
            --bfile ${plink_chrX_basename} \
            --make-bed --out ${plink_chrX_QCed_basename} \
            --freqdiff-x ${bonfx} > XwasQC.log
        """
}