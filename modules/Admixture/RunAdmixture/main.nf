process RunAdmixture {
    publishDir "${params.outdir}/admixture/", mode: 'copy'

    input:
        tuple val(pruned_basename), path(pruned_files)

    output:
        path(admixture_table), emit: admixture_table
        path(admixture_log)
    
    script:
        admixture_table = "${pruned_basename}.${params.admixture_K}.Q"
        admixture_log = "admixture_K${params.admixture_K}.log"
        
        """
        set -eo pipefail

        ${params.tools.admixture} --cv ${pruned_basename}.bed ${params.admixture_K} > ${admixture_log}; 
        """
}
