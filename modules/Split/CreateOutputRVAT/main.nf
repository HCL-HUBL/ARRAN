process CreateOutputRVAT {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "split/$it" }, mode: 'copy'

    label 'GWAS'
    cpus 1
    memory { 100.MB * task.attempt }
    time { 5.minute * task.attempt }
    
    input:
        tuple val(baseqced_basename), path(baseqced_files)
        path(saige_regions)

    output:
        tuple val("${baseqced_basename}_RVAT"), path("${baseqced_basename}_RVAT.{bim,bed,fam}"), emit: plink_RVAT
        path("CreateOutputRVAT.log")

    script:
        extract_cmd = ""
        if(saige_regions) extract_cmd = "--extract range ${params.saige_regions}"

        """
        set -eo pipefail

        ${params.tools.plink} \
            --bfile ${baseqced_basename} \
            --allow-no-sex \
            ${extract_cmd} \
            --make-bed \
            --out ${baseqced_basename}_RVAT > CreateOutputRVAT.log
        """
}