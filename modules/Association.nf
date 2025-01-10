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
        if(params.saige_trait == "binary") binary_flag = "-b TRUE"

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

    output:
        path(gmmat_file), emit: gmmat
        path(vr_file), emit: vr

    script:
        gmmat_file = "${plink_basename}_saige.rda"
        vr_file    = "${plink_basename}_saige.varianceRatio.txt"

        invnorm = "--invNormalize=FALSE"
        if(params.saige_trait == "quantitative") invnorm = "--invNormalize=TRUE"

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
            --traitType=${params.saige_trait} \
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
        path(saige_sv_output), emit: saige_sv_output

    script:
        saige_sv_output = "${plink_basename}_saige.single_variant.tsv"

        """
        set -eo pipefail

        ${params.tools.Rscript} ${params.tools.saige_folder}/step2_SPAtests.R \
            --bedFile=${plink_basename}.bed \
            --bimFile=${plink_basename}.bim \
            --famFile=${plink_basename}.fam \
            --GMMATmodelFile=${gmmat_file} \
            --varianceRatioFile=${vr_file} \
            --is_Firth_beta=TRUE --pCutoffforFirth=0.05 \
            --LOCO=FALSE \
            --is_output_markerList_in_groupTest=TRUE \
            --SAIGEOutputFile=${saige_sv_output}
        """
}
