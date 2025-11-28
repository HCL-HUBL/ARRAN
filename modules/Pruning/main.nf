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