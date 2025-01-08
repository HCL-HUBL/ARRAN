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

        """
        set -eo pipefail

        ${params.tools.Rscript} ${projectDir}/bin/create_phenoFile.R \
            -f ${plink_basename}.fam \
            -c ${covar_file}
        """
}

process GWASFitNullModel {
    input:
        tuple val(plink_basename), path(plink_files)
        path(phenofile)

    output:
        path(GMMATmodel)

    script:
        GMMATmodel = "${plink_basename}_saige.rda"

        invnorm = "--invNormalize=FALSE"
        if(params.saige_trait == "quantitative") invnorm = "--invNormalize=TRUE"

        """
        set -eo pipefail

        ${params.tools.Rscript} ${params.tools.step1_fitNULLGLMM} \
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