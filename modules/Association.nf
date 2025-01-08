#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process GWASFitNullModel {
    input:
        tuple val(plink_basename), path(plink_files)
        path(phenofile)

    output:
        path(GMMATmodel)

    script:
        GMMATmodel = "${plink_basename}_saige.rda"
        params.trait == "quantitative" ? invnorm = "--invNormalize=TRUE" : invnorm = "--invNormalize=FALSE"

        """
        set -eo pipefail

        ${params.tools.Rscript} ${params.tools.SAIGE}/extdata/step1_fitNULLGLMM.R \
            --plinkFile ${plink_basename} \
            --phenoFile ${phenofile} \
            --phenoCol=PHENOTYPE \
            --covarColList=${params.saige_covar} \
            --qCovarColList=${params.saige_qcovar} \
            --sampleIDColinphenoFile=IID \
            ${invnorm} \
            --traitType=${params.rvat_trait} \
            --outputPrefix=${plink_basename}_saige \
            --nThreads=${task.cpus} \
            --IsOverwriteVarianceRatioFile=TRUE
        """
}