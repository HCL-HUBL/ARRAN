# ARRAN

ARRAN (Automatic and Reproducible Rare variants and gwas Analyses in Nextflow) is a Nextflow pipeline designed to perform Genome Wide Association Studies (GWAS) and Rare Variants Association Tests (RVAT).

This pipeline uses Plink to perform QC steps and SAIGE/XWAS to perform the association tests. It includes:
 
 - GWAS and Rare Variants Association Tests for Binary or Continuous traits

 - Inclusion of chrX (QC specific steps and stratified association tests)

 - Generation of logs and plots for the different steps of the pipeline (PCA, Manhattan, QQplot...)

Please refer to our [wiki](https://github.com/HCL-HUBL/ARRAN/wiki) for detailed information about installing, running and interpreting the output of ARRAN.nf.
