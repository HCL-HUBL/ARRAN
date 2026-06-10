process RegionPlot {
    publishDir "${params.outdir}/plots/", mode: 'copy'

    label 'GWAS'
    cpus 1
    memory { 1.GB * task.attempt }
    time { 5.minute * task.attempt }
    
    input:
        tuple val(plink_basename), path(plink_files)
        path(summ_stats)
        path(glist)

    output:
        path("regionplot.pdf")

    script:
        """
        set -eo pipefail

        # To get the rsid of the smallest p-value:
        rsid=\$(awk 'NR == 1 || \$8<min {rsid=\$9; min=\$8} END { print rsid }' ${summ_stats});

        ${params.tools.plink} \
            --bfile ${plink_basename} \
            --r2 --ld-snp \${rsid} --ld-window-kb 500 --ld-window 200 --ld-window-r2 0 \
            --out \${rsid}
        # Converting all consecutive spaces into one tab:
        cat \${rsid}.ld | tr -s ' ' '\t' > \${rsid}.ld.tsv

        regionplot.R \
            -i ${summ_stats} \
            -s \${rsid} \
            -w 250000 \
            -l \${rsid}.ld.tsv \
            -g ${glist} \
            -o regionplot \
            -v FALSE
        """
}