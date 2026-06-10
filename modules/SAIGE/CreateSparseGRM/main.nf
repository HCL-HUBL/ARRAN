// Creates the SparseGRM matrix for information (is not used in subsequent steps)
// useful to check the kinships within the cohort.
// Needs at least 2000 variants with a MAF > 0.01
process CreateSparseGRM {
    publishDir "${params.outdir}/saige/", mode: 'copy'

    label 'GWAS'
    cpus 10
    memory { 1.GB * task.attempt }
    time { 5.minute * task.attempt }

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