process CreateOutputBaseQC {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "QC/$it" }, mode: 'copy'

    label 'QC'
    cpus 1
    memory { 500.MB * task.attempt }
    time { 5.minute * task.attempt }

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