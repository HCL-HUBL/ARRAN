#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CreatePhenoFile {
    input:
        tuple val(plink_basename), path(plink_files)
        path(covar_file)
    
    output:
        path(phenoFile), emit: phenoFile
    
    script:
        phenoFile = "saige_phenofile.tsv"

        binary_flag = "-b FALSE"
        if(params.trait_type == "binary") binary_flag = "-b TRUE"

        """
        set -eo pipefail

        ${params.tools.Rscript} ${projectDir}/bin/create_phenoFile.R \
            -f ${plink_basename}.fam \
            ${binary_flag} \
            -c ${covar_file}
        """
}

process CreateSparseGRM {
    publishDir "${params.outdir}/saige/", mode: 'copy'
    
    input:
        tuple val(plink_pruned_basename), path(plink_pruned_files)

    output:
        path(sampleIDs), emit: sampleIDs
        path(sparseGRM), emit: sparseGRM

    script:
        sampleIDs = "${plink_pruned_basename}*.sparseGRM.mtx.sampleIDs.txt"
        sparseGRM = "${plink_pruned_basename}*.sparseGRM.mtx"

        """
        set -eo pipefail

        ${params.tools.Rscript} ${params.tools.saige_folder}/createSparseGRM.R       \
            --plinkFile=${plink_pruned_basename} \
            --nThreads=${task.cpus} \
            --outputPrefix=${plink_pruned_basename} \
            --numRandomMarkerforSparseKin=2000 \
            --relatednessCutoff=0.125
        """
}

process SaigeFitNullModel {
    publishDir "${params.outdir}/saige/", mode: 'copy'

    input:
        tuple val(plink_basename), path(plink_files)
        path(phenofile)
        val(step)  // should be "GWAS" or "RVAT" will determine if 'isCateVarianceRatio' is TRUE or FALSE

    output:
        path(gmmat_file), emit: gmmat
        path(vr_file), emit: vr

    script:
        gmmat_file = "${plink_basename}_saige.rda"
        vr_file    = "${plink_basename}_saige.varianceRatio.txt"

        invnorm = "--invNormalize=FALSE"
        if(params.trait_type == "quantitative") invnorm = "--invNormalize=TRUE"

        cateVR_cmd = ""
        if(step == "GWAS") cateVR_cmd = "--isCateVarianceRatio=FALSE"
        if(step == "RVAT") cateVR_cmd = "--isCateVarianceRatio=TRUE"

        """
        set -eo pipefail

        ${params.tools.Rscript} ${params.tools.saige_folder}/step1_fitNULLGLMM.R \
            --plinkFile ${plink_basename} \
            --phenoFile ${phenofile} \
            --phenoCol=PHENOTYPE \
            --covarColList=${params.saige_covar} \
            --qCovarColList=${params.saige_qcovar} \
            --sampleIDColinphenoFile=IID \
            ${invnorm} \
            ${cateVR_cmd} \
            --traitType=${params.trait_type} \
            --outputPrefix=${plink_basename}_saige \
            --nThreads=${task.cpus} \
            --IsOverwriteVarianceRatioFile=TRUE
        """
}

process SaigeSingleAssoc {
    publishDir "${params.outdir}/saige/", mode: 'copy'

    input:
        tuple val(plink_basename), path(plink_files)
        path(gmmat_file)
        path(vr_file)
    
    output:
        path(saige_sv), emit: saige_sv

    script:
        saige_sv = "${plink_basename}_saige.single_variant.tsv"

        xpar_cmd = "--is_rewrite_XnonPAR_forMales=FALSE"
        if(params.genome_build == "hg19") xpar_cmd = "--X_PARregion='60001-2699520,154931044-155260560' --is_rewrite_XnonPAR_forMales=TRUE --sampleFile_male=males.list"
        if(params.genome_build == "hg38") xpar_cmd = "--X_PARregion='10001-2781479,155701383-156030895' --is_rewrite_XnonPAR_forMales=TRUE --sampleFile_male=males.list"

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
            ${xpar_cmd} \
            --LOCO=FALSE \
            --is_output_markerList_in_groupTest=TRUE \
            --SAIGEOutputFile=${saige_sv}
        """
}

process ManhattanPlot {
    publishDir "${params.outdir}/plots/", mode: 'copy'

    input:
        path(saige_sv)

    output:
        path(manhattan)
    
    script:
        manhattan = "Manhattan_${saige_sv.baseName}.pdf"
        
        """
        set -eo pipefail

        ${params.tools.Rscript} ${projectDir}/bin/manhattan_plot.R -i ${saige_sv} -o ${manhattan}
        """
}
