#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Processes to split the data into autosomes and chrX and to split each subset 
// into files ready for single-variants analyses and gene-based analyses.

// Process that will split the SNPs belonging to the autosomes (1-22 + 25(PAR)) vs. chrX (23)
process Split_Autosomes_ChrX {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "split/$it" }, mode: 'copy'

    input:
        tuple val(baseqced_basename), path(baseqced_files)

    output:
        tuple val(autosomes_basename), path(autosomes_files), emit: autosomes, optional: true
        tuple val(x_basename), path(x_files), emit: chrx, optional: true

    script:
        autosomes_basename = "${params.out_basename}_autosomes"
        autosomes_files    = "${params.out_basename}_autosomes.{bim,bed,fam}"

        x_basename = "${params.out_basename}_chrX"
        x_files    = "${params.out_basename}_chrX.{bim,bed,fam}"

        """
        set -eo pipefail
        
        awk '{if((\$1 >= 1 && \$1 <= 22) || (\$1 == 25)) { print \$2 }}' ${baseqced_basename}.bim > autosomes_tmp;
        awk '{if(\$1 == 23) { print \$2 }}' ${baseqced_basename}.bim > chrX_tmp;

        if [ `cat autosomes_tmp | wc -l` -gt 0 ]; then
            ${params.tools.plink} \
                --bfile ${baseqced_basename} \
                --allow-no-sex \
                --chr 1-22, 25 \
                --make-bed \
                --out ${autosomes_basename}
        fi

        if [ `cat chrX_tmp | wc -l` -gt 0 ]; then
            ${params.tools.plink} \
                --bfile ${baseqced_basename} \
                --allow-no-sex \
                --chr 23 \
                --make-bed \
                --out ${x_basename}
        fi
        """
}


process CreateOutputGWAS {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "split/$it" }, mode: 'copy'

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


process CreateOutputRVAT {
    publishDir "${params.outdir}/", saveAs: { it.endsWith(".log") ? "logs/$it" : "split/$it" }, mode: 'copy'

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
            --max-maf ${params.rvat_maf} \
            --allow-no-sex \
            ${extract_cmd} \
            --make-bed \
            --out ${baseqced_basename}_RVAT > CreateOutputRVAT.log
        """
}
