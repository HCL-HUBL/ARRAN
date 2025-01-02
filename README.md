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

## The Pipeline

Setting up the environment:

```shell
PLINK="/srv/scratch/chu-lyon.fr/molitorco/scripts/plink1.9/plink"

input_folder="./test_dataset/"
input_basename="EUR_height"
out_folder="./test_dataset/results/"
```

### Input data

To run this pipeline you will need genotyping data in the *plink* format:
 
 - [.bed](https://www.cog-genomics.org/plink/1.9/formats#bed): file representing the genotype calls. Must be accompanied by .bim and .fam formats.

 - [.bim](https://www.cog-genomics.org/plink/1.9/formats#bim): file listing the variants (position and alleles).
  
 - [.fam](https://www.cog-genomics.org/plink/1.9/formats#fam): file with the samples information (ID, family IDs, sex and phenotype).

More information about the formats can be found in the Plink documentation: https://www.cog-genomics.org/plink/1.9/formats.

A test dataset of simulated height in European samples from the 1000 Genomes Project is available in "./test_dataset/".

### GWAS Quality Control

#### Standard QC on the genotype data:

The first step is to filter individuals and variants with low genotyping rates:

 - If needed: remove indivuals that must be filtered, for example because of bad quality or removed consent (--remove)

 - Remove individuals with 5% missing genotypes (--mind 0.05)

 - Remove variants with 5% missing genotypes (--geno 0.05)

 - Produce the eigenvectors and eigenvalues of the genomic PCA

This will generated a new set of plink files with the *_QCed* suffix (--out).

```shell 
${PLINK} \
    --bfile ${input_folder}/${input_basename} \
    --geno 0.05 \
    --mind 0.05 \
    --pca \
    --allow-no-sex \
    --make-bed \
    --out ${out_folder}/${input_basename}_QCed
    # --remove /path/to/samples_to_exclude.txt \
```

Plink will output a log detailing the filtering, in the case of the test dataset:
 
 - 14 people removed due to missing genotype data (--mind).

 - 1 variant removed due to missing genotype data (--geno).

 - Total genotyping rate in remaining samples is 0.999816.

The PCA results are written to: 

 - ./test_dataset/results/EUR_height_QCed.eigenval

 - ./test_dataset/results/EUR_height_QCed.eigenvec

#### Looking at the PCA:

```shell
Rscript ./bin/PCA.R -i ${out_folder}/${input_basename}_QCed.eigenvec -f ${out_folder}/${input_basename}_QCed.fam -v
```


### GWAS association 

### Downstream analysis

