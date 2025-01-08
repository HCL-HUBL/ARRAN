# HCL-GWAS

## Introduction

Nextflow pipeline to perform Genome Wide Association Studies (GWAS) and/or Rare Variants Association Tests (RVAT). This Readme contains step-by-step instructions to configure and run the pipeline.

The pipeline works for binary and continuous traits and uses SAIGE to perform the association tests.

## Table of Contents

1. [Introduction](#introduction)
2. [Table of Contents](#table-of-contents)
3. [Dependencies](#dependencies)
3. [The Pipeline](#the-pipeline)
    - [Input data](#input-data)
    - [Quality Control](#quality-control)
    - [GWAS association](#gwas-association)
    - [Downstream analysis](#downstream-analysis)

## Dependencies

The following tools need to be installed on your machine: 

 - [plink v1.9](https://www.cog-genomics.org/plink/1.9/)

 - R 
    - The R packages 'optparse' and 'ggplot2'

 - SAIGE 

All dependencies are available inside a Singularity image, that can be build from the recipe provided within this repository [HCL-GWAS.def](./HCL-GWAS.def):

```shell
singularity build HCL-GWAS.sif HCL-GWAS.def
```

*Note: The scripts were developed and tested on Linux (Debian release 11) using nextflow v22.04.5 and R v4.4.1*

## The Pipeline

The pipeline can be launched from [HCL-GWAS.nf](./HCL-GWAS.nf).

You will need to change values in the configuration file [default.conf](./confs/default.conf) to adjust the QC and association steps to your own needs.

### Input data

To run this pipeline you will need genotyping data in the *plink* format:
 
 - [.bed](https://www.cog-genomics.org/plink/1.9/formats#bed): file representing the genotype calls. Must be accompanied by .bim and .fam formats.

 - [.bim](https://www.cog-genomics.org/plink/1.9/formats#bim): file listing the variants (position and alleles).
  
 - [.fam](https://www.cog-genomics.org/plink/1.9/formats#fam): file with the samples information (ID, family IDs, sex and phenotype).

More information about the formats can be found in the Plink documentation: https://www.cog-genomics.org/plink/1.9/formats.

You will also need a [covariates file](https://www.cog-genomics.org/plink/1.9/input#covar) which **must** have a header.

### Quality Control

#### QC on the genotype data:

The first step filters performs standard GWAS quality control:

 - Remove individuals with >5% missing genotypes (can be changed with 'qc_mind')

 - Remove individuals listed in the file designed by 'qc_remove' (keep *qc_remove = ""* if you do not wish to filter any individual)

 - Remove individuals with extreme heterozygosity ('qc_hetfilter'):
    - 'low':  remove individuals with a F coefficient > 3 SDs from the cohort's mean.
    - 'high': remove individuals with a F coefficient < 3 SDs from the cohort's mean.
    - 'both': applies both of the above.
    - *note*: the F coefficient is inversely correlated with heterozygosity (so a 'high' heterozygosity corresponds to a low F value).

 - Remove variants with 5% missing genotypes (can be changed with 'qc_geno')

 - Produce the eigenvectors and eigenvalues of the genomic PCA

 - Remove variants based on their Minor Allele Frequency:
    - For the GWAS: variants with a MAF < 'qc_maf'
    - For the RVAT: variants with a MAF > 'rvat_maf'

This pipeline will output the following files:

 - <basename>_QCed.{bim,bed,fam}: contains the variants and individuals that passed the QC step

 - ./plots/

### GWAS association 

### Downstream analysis

