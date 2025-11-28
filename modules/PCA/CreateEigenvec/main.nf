process CreateEigenvec {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "QC/$it" }, mode: 'copy'

    input:
        tuple val(baseqced_basename), path(baseqced_files)
    
    output:
        path("${eigen_basename}.eigenval"), emit: eigenval
        path("${eigen_basename}.eigenvec"), emit: eigenvec
        path("CreateEigenvec.log")
    
    script:
        eigen_basename = "${baseqced_basename}_PCA"
        eigen_files    = "${baseqced_basename}_PCA.{bim,bed,fam}"

        """
        set -eo pipefail

        ${params.tools.plink} \
            --bfile ${baseqced_basename} \
            --maf ${params.gwas_maf} \
            --allow-no-sex \
            --pca \
            --make-bed \
            --out ${eigen_basename} > CreateEigenvec.log
        """
}