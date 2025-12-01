process HetFilter {
    publishDir "${params.outdir}/", saveAs: { it.endsWith("valides") ? "QC/$it" : "plots/$it" }, mode: 'copy'

    input:
        path(geno_het)
    
    output:
        path("${geno_het}.valides"), emit: valides
        path("${geno_het}.nonvalides"), emit: nonvalides
        path("${geno_het}.pdf")

    script:
        """
        set -eo pipefail
        
        het_check.R -i ${geno_het} -f ${params.qc_hetfilter}
        """
}
