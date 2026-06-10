// Perform the Single Variant Association Test (= GWAS) using SAIGE+
process SaigeSingleAssoc {
    publishDir "${params.outdir}/saige/", mode: 'copy'

    label 'GWAS'
    cpus 1
    memory { 1.GB * task.attempt }
    time { 1.hour * task.attempt }

    input:
        tuple val(plink_basename), path(plink_files)
        path(gmmat_file)
        path(vr_file)
    
    output:
        path(saige_sv), emit: saige_sv

    script:
        saige_sv = "${plink_basename}_saige.single_variant.tsv"

        // xpar_cmd = "--is_rewrite_XnonPAR_forMales=FALSE"
        // if(params.genome_build == "hg19") xpar_cmd = "--X_PARregion='60001-2699520,154931044-155260560' --is_rewrite_XnonPAR_forMales=TRUE --sampleFile_male=males.list"
        // if(params.genome_build == "hg38") xpar_cmd = "--X_PARregion='10001-2781479,155701383-156030895' --is_rewrite_XnonPAR_forMales=TRUE --sampleFile_male=males.list"

        """
        set -eo pipefail

        awk '\$5 == 1 { print \$2 }' ${plink_basename}.fam > males.list

        ${params.tools.Rscript} ${params.tools.saige_folder}/step2_SPAtests.R \
            --bedFile=${plink_basename}.bed \
            --bimFile=${plink_basename}.bim \
            --famFile=${plink_basename}.fam \
            --GMMATmodelFile=${gmmat_file} \
            --varianceRatioFile=${vr_file} \
            --is_Firth_beta=TRUE --pCutoffforFirth=0.05 \
            --LOCO=FALSE \
            --is_output_markerList_in_groupTest=TRUE \
            --SAIGEOutputFile=${saige_sv}
        """
}