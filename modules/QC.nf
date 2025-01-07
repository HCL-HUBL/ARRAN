#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process BaseQC {
    input:
        tuple val(plink_basename), path(plink_files)
        path(qc_remove)

    output:
        tuple val("${plink_basename}_baseQC"), path("${plink_basename}_baseQC.bim"), path("${plink_basename}_baseQC.bed"), path("${plink_basename}_baseQC.fam"), emit: plink_baseQC
        path("${plink_basename}_baseQC.eigenvec"), emit: eigenvec
        path("${plink_basename}_baseQC.eigenval"), emit: eigenval
        path("BaseQC.log")

    script:
        remove_cmd = ""
        if(qc_remove) remove_cmd = "--remove ${params.qc_remove}"

        """
        set -eo pipefail

        ${params.tools.plink} \
            --bfile ${plink_basename} \
            --geno ${params.qc_geno} \
            --mind ${params.qc_mind} \
            --maf  ${params.qc_maf} \
            ${remove_cmd} \
            --pca \
            --allow-no-sex \
            --make-bed \
            --out ${plink_basename}_baseQC > BaseQC.log
        """
}

process PlotPCA {
    publishDir "${params.outdir}/plots/", mode: 'copy'

    input:
        tuple val(baseqc_basename), path("${baseqc_basename}.bim"), path("${baseqc_basename}.bed"), path("${baseqc_basename}.fam")
        path(eigenvec)

    output:
        path("${baseqc_basename}_PCA.pdf")

    script:
        """
        set -eo pipefail

        ${params.tools.Rscript} ${projectDir}/bin/PCA.R \
            -i ${eigenvec} \
            -f ${baseqc_basename}.fam \
            -o ${baseqc_basename}_PCA
        """
}

process Pruning {
    publishDir "${params.outdir}/", mode: 'copy'

    input:
        tuple val(baseqc_basename), path("${baseqc_basename}.bim"), path("${baseqc_basename}.bed"), path("${baseqc_basename}.fam")

    output:
        path("${baseqc_basename}.prune.in"), emit: prune_in
        path("${baseqc_basename}.prune.out"), emit: prune_out
        path("Pruning.log")

    script:
        """
        set -eo pipefail

        ${params.tools.plink} \
            --bfile ${baseqc_basename} \
            --indep-pairwise ${params.pr_window} ${params.pr_step} ${params.pr_r2} \
            --out ${baseqc_basename} > Pruning.log
        """
}

process HetCoeff {
    publishDir "${params.outdir}/", mode: 'copy'

    input:
        tuple val(baseqc_basename), path("${baseqc_basename}.bim"), path("${baseqc_basename}.bed"), path("${baseqc_basename}.fam")
        path(prune_in)

    output:
        path("${baseqc_basename}.het"), emit: het
        path("HetCoeff.log")

    script:
        """
        set -eo pipefail
        
        ${params.tools.plink} \
            --bfile ${baseqc_basename} \
            --extract ${prune_in} \
            --het \
            --out ${baseqc_basename} > HetCoeff.log
        """
}

process HetFilter {
    publishDir "${params.outdir}/", saveAs: { it.endsWith("valides") ? "$it" : "plots/$it" }, mode: 'copy'

    input:
        path(geno_het)
    
    output:
        path("${geno_het}.valides"), emit: valides
        path("${geno_het}.nonvalides"), emit: nonvalides
        path("${geno_het}.pdf")

    script:
        """
        set -eo pipefail

        ${params.tools.Rscript} ${projectDir}/bin/het_check.R \
            -i ${geno_het} \
            -f ${params.qc_hetfilter}
        """
}

process CreateOutput {
    publishDir "${params.outdir}/", mode: 'copy'

    input:
        tuple val(baseqc_basename), path("${baseqc_basename}.bim"), path("${baseqc_basename}.bed"), path("${baseqc_basename}.fam")
        path(valides)

    output:
        tuple val("${baseqc_basename}ed"), path("${baseqc_basename}ed.bim"), path("${baseqc_basename}ed.bed"), path("${baseqc_basename}ed.fam"), emit: plink_QCed
        
    script:
        """
        set -eo pipefail

        ${params.tools.plink} \
            --bfile ${baseqc_basename} \
            --keep ${valides} \
            --make-bed \
            --out ${baseqc_basename}ed
        """
}