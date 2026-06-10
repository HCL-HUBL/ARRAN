process ManhattanPlot {
    publishDir "${params.outdir}/plots/", mode: 'copy'

    label 'GWAS'
    cpus 1
    memory { 1.GB * task.attempt }
    time { 5.minute * task.attempt }

    input:
        path(assoc_tsv)
        val(chr_col)
        val(pos_col)
        val(marker_col)
        val(pval_col)

    output:
        path(manhattan), optional: true
    
    script:
        manhattan = "Manhattan_${assoc_tsv.baseName}_${pval_col}.pdf"
        
        """
        set -eo pipefail

        manhattan_plot.R \
            -i ${assoc_tsv} \
            -c ${chr_col} -b ${pos_col} -m ${marker_col} -p ${pval_col} \
            -o ${manhattan}
        """
}