// Creates the Group File for the Rare Variant Association Test with SAIGE+
// This contains the list of genes (or regions) to consider, with their associated variants.
// We are using glist-hg19 or glist-hg38 to get the genes boundaries and assign variants.
process CreateGroupFile {
    publishDir "${params.outdir}/saige/", mode: 'copy'

    label 'GWAS'
    cpus 1
    memory { 2.GB * task.attempt }
    time { 30.minute * task.attempt }

    input:
        tuple val(plink_basename), path(plink_files)
        path(glist)
        path(annot)
    
    output:
        path(groupFile), emit: group_file

    script:
        groupFile = "saige.groupFile"

        annot_cmd = ""
        if(annot) annot_cmd = "-a ${annot}"

        """
        set -eo pipefail

        echo -e "CHR\\tSNP\\tBP" > variants.report
        cut -f 1,2,4 ${plink_basename}.bim >> variants.report
        
        ${params.tools.plink} --annotate variants.report ranges=${glist} --out genes

        awk -F' ' '{split(\$4, genes, "|"); for (i in genes) {print \$1, \$2, \$3, genes[i]}}' genes.annot > genes.annot.split
        sed 's/\t/ /g' genes.annot.split | cut -d " " -f 2,4,5 | sed 's/(0)//g' > genes.annot.clean

        create_groupFile.R -g genes.annot.clean ${annot_cmd}
        """
}