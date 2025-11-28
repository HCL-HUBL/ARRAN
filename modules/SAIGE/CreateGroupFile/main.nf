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