process BaseQC {
    publishDir "${params.outdir}/logs/",  pattern: "BaseQC.log", mode: 'copy'

    label 'QC'
    cpus 1
    memory { 5.GB * task.attempt }
    time { 5.minute * task.attempt }

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