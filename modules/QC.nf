#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// This process will split the PAR regions into an separate chr if it does not exists
// NON-PAR genotypes on chrX for males will be doubled
// Unknown variant IDs "." will be transformed to "chr_pos_ref_alt"
process GenotypesPreprocessing {
    input:
        tuple val(plink_basename), path(plink_files)

    output:
        tuple val(preprocessed_basename), path(preprocessed_files), emit: plink_preprocessed

    script:
        preprocessed_basename = "${params.out_basename}_base"
        preprocessed_files    = "${params.out_basename}_base.{bim,bed,fam}"

        """
        set -eo pipefail

        # Replacing missing variant IDs '.' with 'chr_pos_ref_alt':
        awk '{if(\$2 == ".") {print  \$1"\\t"\$1"_"\$4"_"\$5"_"\$6"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6} else {print \$0}}' ${plink_basename}.bim > tmp && mv tmp ${plink_basename}.bim

        # If chrXY (PAR) is not present, we create it:
        if ! grep -q "^25" ${plink_basename}.bim;
        then
            ${params.tools.plink} \
                --bfile ${plink_basename} \
                --split-x '${params.genome_build}' 'no-fail' \
                --allow-no-sex \
                --make-bed \
                --out ${preprocessed_basename}
        else
            ${params.tools.plink} \
                --bfile ${plink_basename} \
                --allow-no-sex \
                --make-bed \
                --out ${preprocessed_basename}
        fi 
        """
}


process BaseQC {
    publishDir "${params.outdir}/logs/",  pattern: "BaseQC.log", mode: 'copy'

    input:
        tuple val(plink_basename), path(plink_files)
        path(qc_remove)

    output:
        tuple val(baseqc_basename), path(baseqc_files), emit: plink_baseQC
        path("BaseQC.log")

    script:
        baseqc_basename = "${plink_basename}QC"
        baseqc_files    = "${plink_basename}QC.{bim,bed,fam}"

        remove_cmd = ""
        if(qc_remove) remove_cmd = "--remove ${params.qc_remove}"

        """
        set -eo pipefail
        
        # Running "geno" first, to avoid removing individuals when there are sets 
        # of variants with high missing % in a subset of samples:
        ${params.tools.plink} \
            --bfile ${plink_basename} \
            ${remove_cmd} \
            --geno ${params.qc_geno} \
            --allow-no-sex \
            --write-snplist \
            --out valid_snps

        # We only perform "mind" on the list of variants extracted before:
        ${params.tools.plink} \
            --bfile ${plink_basename} \
           	--extract valid_snps.snplist \
            ${remove_cmd} \
            --mind ${params.qc_mind} \
            --allow-no-sex \
            --make-bed \
            --out ${baseqc_basename} > BaseQC.log
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


process RunAdmixture {
    publishDir "${params.outdir}/admixture/", mode: 'copy'

    input:
        tuple val(pruned_basename), path(pruned_files)

    output:
        path(admixture_table), emit: admixture_table
        path(admixture_log)
    
    script:
        admixture_table = "${pruned_basename}.${params.admixture_K}.Q"
        admixture_log = "admixture_K${params.admixture_K}.log"
        
        """
        set -eo pipefail

        ${params.tools.admixture} --cv ${pruned_basename}.bed ${params.admixture_K} > ${admixture_log}; 
        """
}


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

        ${params.tools.Rscript} ${projectDir}/bin/admixture.R \
            -i ${admixture_table} \
            -k ${params.admixture_K} \
            -o ${admixture_barplot}
        """
}


process CreateOutputBaseQC {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "QC/$it" }, mode: 'copy'

    input:
        tuple val(baseqc_basename), path(baseqc_files)
        path(prune_in)
        path(valides)

    output:
        tuple val(qced_basename), path(qced_files), emit: plink_QCed
        tuple val(pruned_basename), path(pruned_files), emit: plink_QCed_pruned
        path("CreateOutputBaseQC.log")

    script:
        qced_basename = "${baseqc_basename}ed"
        qced_files    = "${baseqc_basename}ed.{bim,bed,fam}"

        pruned_basename = "${baseqc_basename}ed_pruned"
        pruned_files    = "${baseqc_basename}ed_pruned.{bim,bed,fam}"

        """
        set -eo pipefail

        # Creating the pruned set:
        ${params.tools.plink} \
            --bfile ${baseqc_basename} \
            --keep ${valides} \
            --extract ${prune_in} \
            --allow-no-sex \
            --make-bed \
            --out ${pruned_basename}

        ${params.tools.plink} \
            --bfile ${baseqc_basename} \
            --keep ${valides} \
            --allow-no-sex \
            --make-bed \
            --out ${qced_basename} > CreateOutputBaseQC.log
        """
}


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


// Process used to flag problematic variants in term of HWE
process HWEFlag {
    errorStrategy 'ignore'

    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "QC/$it" }, mode: 'copy'

    input:
        tuple val(baseqc_basename), path(baseqc_files)
        val(qc_hwe)

    output:
        path(hwelist), optional: true
        path("HWEFlag.log")
    
    script:
        hwelist = "hwe_filtered.snplist"

        """
        set -eo pipefail

        ${params.tools.plink} \
            --allow-no-sex \
            --bfile ${baseqc_basename} \
            --hwe ${qc_hwe} \
            --write-snplist \
            --out hwe_pass > HWEFlag.log

        # We exclude variants which passed hwe threshold, 
        # to only keep those which didn't pass the filter:
        ${params.tools.plink} \
            --allow-no-sex \
            --bfile ${baseqc_basename} \
            --exclude hwe_pass.snplist \
            --write-snplist \
            --out hwe_filtered
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
        binary_cmd = "-b FALSE"
        if(params.trait_type == "binary") binary_cmd = "-b TRUE"

        """
        set -eo pipefail

        ${params.tools.Rscript} ${projectDir}/bin/PCA.R \
            -i ${eigenvec} \
            -f ${qced_basename}.fam \
            ${binary_cmd} \
            -o ${qced_basename}_PCA
        """
}
