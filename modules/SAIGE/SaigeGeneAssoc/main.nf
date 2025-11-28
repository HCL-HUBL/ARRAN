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