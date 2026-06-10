// Fits the Null Model for the GWAS and/or RVAT. This model will be compared with 
// the full model (containing the genotypes) for association testing.
process SaigeFitNullModel {
    publishDir "${params.outdir}/saige/", mode: 'copy'

    label 'GWAS'
    cpus 10
    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }

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
            --isCovariateOffset=FALSE \
            ${invnorm} \
            ${cateVR_cmd} \
            --traitType=${params.trait_type} \
            --outputPrefix=${plink_basename}_saige \
            --nThreads=${task.cpus} \
            --IsOverwriteVarianceRatioFile=TRUE
        """
}