process CNVFiltering {
    publishDir "${params.outdir}/cnvs", mode: 'copy'

    input:
        path(cnv_file)

    output:
        path "${cnv_file}.filtered"

    script:
        sample = cnv_file.baseName
        
        """
        # Filtrer les CNVs par taille (e.g., supprimer ceux < 1kb)
        awk '\$6 >= 1000' ${cnv_file} > ${sample}_filtered.cnv

        # Annoter les CNVs avec les gènes (exemple simplifié)
        # Supposons que gene_annotation_file est un fichier BED avec les coordonnées des gènes
        bedtools intersect -a ${sample}_filtered.cnv -b ${gene_annotation_file} -wa -wb > ${sample}_filtered_annotated.cnv

        # Nettoyer les fichiers intermédiaires
        rm ${sample}_filtered.cnv
        """
}