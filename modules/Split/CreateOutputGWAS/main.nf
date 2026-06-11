process CreateOutputGWAS {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "split/$it" }, mode: 'copy'

    label 'GWAS'
    cpus 1
    memory { 500.MB * task.attempt }
    time { 5.minute * task.attempt }

    input:
        tuple val(baseqced_basename), path(baseqced_files)
        path(saige_regions)

    output:
        tuple val(gwas_basename), path(gwas_files), emit: plink_GWAS
        path("CreateOutputGWAS.log")

    script:
        gwas_basename = "${baseqced_basename}_GWAS"
        gwas_files    = "${baseqced_basename}_GWAS.{bim,bed,fam}"

        extract_cmd = ""
        if(saige_regions) extract_cmd = "--extract range ${params.saige_regions}"

        """
        set -eo pipefail

        ${params.tools.plink} \
            --bfile ${baseqced_basename} \
            --maf ${params.gwas_maf} \
            --allow-no-sex \
            ${extract_cmd} \
            --make-bed \
            --out ${gwas_basename} > CreateOutputGWAS.log
        """
}
