// Creates the phenotype file for SAIGE which contains the samples IIDs and 
// the phenotypes + covariates 
process CreatePhenoFile {
    publishDir "${params.outdir}/saige/", mode: 'copy'

    label 'GWAS'
    cpus 1
    memory { 100.MB * task.attempt }
    time { 5.minute * task.attempt }

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

        create_phenoFile.R \
            -f ${plink_basename}.fam \
            -e ${eigenvec} \
            -c ${covar_file} \
            ${binary_flag}
        """
}