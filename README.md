# HCL-GWAS

## Introduction

This readme contains step-by-step instructions to perform a Genome Wide Association Study (GWAS):

- [QC Steps](#gwas-quality-control)

- [GWAS Association](#gwas-association)

- [Downstream Analysis](#downstream-analysis)

## Table of Contents

1. [Introduction](#introduction)
2. [Table of Contents](#table-of-contents)
3. [Dependencies](#dependencies)
3. [The Pipeline](#the-pipeline)
    - [Input data](#input-data)
    - [GWAS Quality Control](#gwas-quality-control)
    - [GWAS association](#gwas-association)
    - [Downstream analysis](#downstream-analysis)

## Dependencies

 - plink v1.9

 - R
 
 - 

## The Pipeline

Setting up the environment:

```shell
PLINK="/path/to/plink"

input_folder="./"
input_basename="Final_dat_merged_corrected"
out_folder="./"

cd ${out_folder}
```

### Input data

To run this pipeline you will need genotyping data in the *plink* format:
 
    - [.bed](https://www.cog-genomics.org/plink/1.9/formats#bed): file representing the genotype calls. Must be accompanied by .bim and .fam formats.

    - [.bim](https://www.cog-genomics.org/plink/1.9/formats#bim): file listing the variants (position and alleles).
  
    - [.fam](https://www.cog-genomics.org/plink/1.9/formats#fam): file with the samples information (ID, family IDs, sex and phenotype).

More information about the formats can be found in the Plink documentation: https://www.cog-genomics.org/plink/1.9/formats.

### GWAS Quality Control

The first step is to filter individuals and variants with low genotyping rates:

    - Remove individuals with 5% missing genotypes (--mind 0.05)

    - Remove variants with 5% missing genotypes (--geno 0.05)

    - Remove any indivual that must be filtered, for example because of duplicates or removed consent (--remove)



```shell 
${PLINK} \
    --bfile ${input_file} \
    --geno 0.05 \
    --mind 0.05 \
    --remove /path/to/samples_to_exclude.txt \
    --allow-no-sex \
    --make-bed \
    --out ${out_folder}/${input_basename}_QCed
```

### GWAS association 

### Downstream analysis

