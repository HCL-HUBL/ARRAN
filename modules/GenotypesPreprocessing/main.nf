// This process will split the PAR regions into an separate chr if it does not exists
// NON-PAR genotypes on chrX for males will be doubled
// Unknown variant IDs "." will be transformed to "chr_pos_ref_alt"
process GenotypesPreprocessing {    
    label 'QC'
    cpus 1
    memory { 5.GB * task.attempt }
    time { 5.minute * task.attempt }

    input:
        tuple val(plink_basename), path(plink_files)

    output:
        tuple val(preprocessed_basename), path(preprocessed_files), emit: plink_preprocessed

    script:
        preprocessed_basename = "${params.out_basename}_base"
        preprocessed_files    = "${params.out_basename}_base.{bim,bed,fam}"

        """
        set -eo pipefail

        # Replacing missing variant IDs '.' with 'chr_pos_ref_alt':
        awk '{if(\$2 == ".") {print  \$1"\\t"\$1"_"\$4"_"\$5"_"\$6"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6} else {print \$0}}' ${plink_basename}.bim > tmp && mv tmp ${plink_basename}.bim

        # If chrXY (PAR) is not present, we create it:
        if ! grep -q "^25" ${plink_basename}.bim;
        then
            ${params.tools.plink} \
                --bfile ${plink_basename} \
                --split-x '${params.genome_build}' 'no-fail' \
                --allow-no-sex \
                --make-bed \
                --out ${preprocessed_basename}
        else
            ${params.tools.plink} \
                --bfile ${plink_basename} \
                --allow-no-sex \
                --make-bed \
                --out ${preprocessed_basename}
        fi 
        """
}