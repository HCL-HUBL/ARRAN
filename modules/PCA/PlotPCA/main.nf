process PlotPCA {
    publishDir "${params.outdir}/plots/", mode: 'copy'

    label 'QC'
    cpus 1
    memory { 500.MB * task.attempt }
    time { 5.minute * task.attempt }

    input:
        tuple val(qced_basename), path(qced_files)
        path(eigenvec)

    output:
        path("${qced_basename}_PCA.pdf")

    script:
        binary_cmd = "-b FALSE"
        if(params.trait_type == "binary") binary_cmd = "-b TRUE"

        """
        set -eo pipefail

        PCA.R \
            -i ${eigenvec} \
            -f ${qced_basename}.fam \
            ${binary_cmd} \
            -o ${qced_basename}_PCA
        """
}