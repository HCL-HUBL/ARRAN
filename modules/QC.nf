#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process BaseQC {
    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple val(plink_basename), path(plink_files)
        path(qc_remove)

    output:
        tuple val("${plink_basename}_QCed"), path("${plink_basename}_QCed.bim"), path("${plink_basename}_QCed.bed"), path("${plink_basename}_QCed.fam"), emit: plink_QCed
        path("${plink_basename}_QCed.eigenvec"), emit: eigenvec
        path("${plink_basename}_QCed.eigenval"), emit: eigenval
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
            --out ${plink_basename}_QCed > BaseQC.log
        """
}


process Pruning {
    publishDir "${params.outdir}/", mode: 'copy'

    input:
        tuple val(qced_basename), path("${qced_basename}.bim"), path("${qced_basename}.bed"), path("${qced_basename}.fam")

    output:
        path("${qced_basename}.prune.in"), emit: prune_in
        path("${qced_basename}.prune.out"), emit: prune_out
        path("Pruning.log")

    script:
    """
    set -eo pipefail

    ${params.tools.plink} \
        --bfile ${qced_basename} \
        --indep-pairwise ${params.pr_window} ${params.pr_step} ${params.pr_r2} \
        --out ${qced_basename} > Pruning.log
    """
}


process HetCoeff {
    publishDir "${params.outdir}/", mode: 'copy'

    input:
        tuple val(qced_basename), path(qced_files)
        path(prune_in)

    output:
        path("${qced_basename}.het"), emit: het
        path("HetCoeff.log")

    script:
    """
    set -eo pipefail
    
    ${params.tools.plink} \
        --bfile ${qced_basename} \
        --extract ${prune_in} \
        --het \
        --out ${qced_basename} > HetCoeff.log
    """
}


process HetFilter {
    publishDir "${params.outdir}/", mode: 'copy'

    input:
        path(geno_het)
    
    output:
        path("${geno_het}.valides"), emit: valides
        path("${geno_het}.nonvalides"), emit: nonvalides
        path("${geno_het}.png")

    script:
    """
    set -eo pipefail

    ${params.tools.Rscript} \
        -i ${geno_het} \
        -f ${params.qc_hetfilter}
    """
}
