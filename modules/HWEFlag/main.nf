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