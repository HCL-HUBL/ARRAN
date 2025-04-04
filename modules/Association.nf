#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Creates the phenotype file for SAIGE which contains the samples IIDs and 
// the phenotypes + covariates 
process CreatePhenoFile {
    publishDir "${params.outdir}/saige/", mode: 'copy'
    
    input:
        tuple val(plink_basename), path(plink_files)
        path(eigenvec)
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
            -e ${eigenvec} \
            -c ${covar_file} \
            ${binary_flag}
        """
}

// Creates the SparseGRM matrix for information (is not used in subsequent steps)
// useful to check the kinships within the cohort.
// Needs at least 2000 variants with a MAF > 0.01
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

// Fits the Null Model for the GWAS and/or RVAT. This model will be compared with 
// the full model (containing the genotypes) for association testing.
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
            --isCovariateOffset=FALSE \
            ${invnorm} \
            ${cateVR_cmd} \
            --traitType=${params.trait_type} \
            --outputPrefix=${plink_basename}_saige \
            --nThreads=${task.cpus} \
            --IsOverwriteVarianceRatioFile=TRUE
        """
}

// Perform the Single Variant Association Test (= GWAS) using SAIGE+
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

// Creates the Group File for the Rare Variant Association Test with SAIGE+
// This contains the list of genes (or regions) to consider, with their associated variants.
// We are using glist-hg19 or glist-hg38 to get the genes boundaries and assign variants.
process CreateGroupFile {
    publishDir "${params.outdir}/saige/", mode: 'copy'

    input:
        tuple val(plink_basename), path(plink_files)
        path(glist)
    
    output:
        path(groupFile), emit: group_file

    script:
        groupFile = "saige.groupFile"

        """
        set -eo pipefail

        echo -e "CHR\\tID\\tBP" > variants.report
        cut -f 1,2,4 ${plink_basename}.bim  >> variants.report
        
        ${params.tools.plink} --annotate variants.report ranges=${glist} 
        sed 's/\\t/ /g' plink.annot  | cut -d " " -f 2,4  | sed 's/(.*//g' > plink.annot.clean

        ${params.tools.Rscript} ${projectDir}/bin/create_groupFile.R -a plink.annot.clean
        """
}

// Perform the Rare Variant Association Tests.
// Only on variants with MAF < 'rvat_maf' (see conf, usually 0.01)
// Using the group file created with the 'CreateGroupFile' process.
process SaigeGeneAssoc {
    publishDir "${params.outdir}/saige/", mode: 'copy'

    input:
        tuple val(plink_basename), path(plink_files)
        path(gmmat_file)
        path(vr_file)
        path(groupFile)
    
    output:
        path(saige_gene), emit: saige_gene

    script:
        saige_gene = "${plink_basename}_saige_gene_based.tsv"
        
        """
        set -eo pipefail

        ${params.tools.Rscript} ${params.tools.saige_folder}/step2_SPAtests.R \
            --bedFile=${plink_basename}.bed \
            --bimFile=${plink_basename}.bim \
            --famFile=${plink_basename}.fam \
            --GMMATmodelFile=${gmmat_file} \
            --varianceRatioFile=${vr_file} \
            --LOCO=FALSE \
            --is_output_markerList_in_groupTest=TRUE \
            --is_single_in_groupTest=FALSE \
            --SAIGEOutputFile=${plink_basename}_saige_gene_based.tsv \
            --groupFile=${groupFile} \
            --annotation_in_groupTest="no_annot" \
            --maxMAF_in_groupTest=0.001,0.01,0.1
        """
// --annotation_in_groupTest="lof,lof:missense,lof:missense:synonymous" \
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

process QQPlot {
    publishDir "${params.outdir}/plots/", mode: 'copy'

    input:
        path(saige_sv)
        val(pcol)

    output:
        path(qqplot)

    script:
        qqplot = "QQplot_${saige_sv.baseName}.pdf"

        """
        set -eo pipefail

        ${params.tools.Rscript} ${projectDir}/bin/QQplot.R -i ${saige_sv} -p ${pcol} -o ${qqplot}
        """
}
