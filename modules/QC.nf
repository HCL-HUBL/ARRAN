#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process BaseQC {
    publishDir "${params.outdir}/logs/",  pattern: "BaseQC.log", mode: 'copy'

    input:
        tuple val(plink_basename), path(plink_files)
        path(qc_remove)

    output:
        tuple val("${plink_basename}_baseQC"), path("${plink_basename}_baseQC.{bim,bed,fam}"), emit: plink_baseQC
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
            ${remove_cmd} \
            --allow-no-sex \
            --make-bed \
            --out ${plink_basename}_baseQC > BaseQC.log
        """
}

process Pruning {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "QC/$it" }, mode: 'copy'

    input:
        tuple val(baseqc_basename), path(baseqc_files)

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
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "QC/$it" }, mode: 'copy'

    input:
        tuple val(baseqc_basename), path(baseqc_files)
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
    publishDir "${params.outdir}/", saveAs: { it.endsWith("valides") ? "QC/$it" : "plots/$it" }, mode: 'copy'

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

process CreateOutputGWAS {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "QC/$it" }, mode: 'copy'

    input:
        tuple val(baseqc_basename), path(baseqc_files)
        path(prune_in)
        path(valides)

    output:
        tuple val("${baseqc_basename}ed"), path("${baseqc_basename}ed.{bim,bed,fam}"), emit: plink_QCed
        tuple val("${baseqc_basename}ed_pruned"), path("${baseqc_basename}ed_pruned.{bim,bed,fam}"), emit: plink_QCed_pruned
        path("${baseqc_basename}ed.eigenval"), emit: eigenval
        path("${baseqc_basename}ed.eigenvec"), emit: eigenvec
        path("CreateOutputGWAS.log")

    script:
        """
        set -eo pipefail

        # Creating the pruned set:
        ${params.tools.plink} \
            --bfile ${baseqc_basename} \
            --keep ${valides} \
            --extract ${prune_in} \
            --maf ${params.gwas_maf} \
            --allow-no-sex \
            --make-bed \
            --out ${baseqc_basename}ed_pruned

        ${params.tools.plink} \
            --bfile ${baseqc_basename} \
            --keep ${valides} \
            --maf ${params.gwas_maf} \
            --allow-no-sex \
            --pca \
            --make-bed \
            --out ${baseqc_basename}ed > CreateOutputGWAS.log
        """
}


process PlotPCA {
    publishDir "${params.outdir}/plots/", mode: 'copy'

    input:
        tuple val(qced_basename), path(qced_files)
        path(eigenvec)

    output:
        path("${qced_basename}_PCA.pdf")

    script:
        """
        set -eo pipefail

        ${params.tools.Rscript} ${projectDir}/bin/PCA.R \
            -i ${eigenvec} \
            -f ${qced_basename}.fam \
            -o ${qced_basename}_PCA
        """
}
