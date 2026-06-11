process HetCoeff {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "QC/$it" }, mode: 'copy'
    
    label 'QC'
    cpus 1
    memory { 500.MB * task.attempt }
    time { 5.minute * task.attempt }
    
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
