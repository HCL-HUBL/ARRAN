#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// Importing Processes from modules:
include { BaseQC }              from './modules/QC.nf'
include { PlotPCA }             from './modules/QC.nf'
include { Pruning }             from './modules/QC.nf'
include { HetCoeff }            from './modules/QC.nf'
include { HetFilter }           from './modules/QC.nf'
include { CreateOutputGWAS }    from './modules/QC.nf'


// Initialising the options with default values:
params.plink_fileset = ""                 // The path to the plink fileset (/path/to/example.{bim,bed,fam})
params.outdir        = "${launchDir}"     // The output directory where the outputs will be stored
params.genome_build  = "hg19"             // The genome build (not used for now)

params.qc_mind       = 0.05               // Individuals with >5% missing genotypes will be removed
params.qc_geno       = 0.05               // Variants with >5% missing genotypes will be removed
params.qc_remove     = ""                 // (optional) File listing IIDs and FIDs of individuals to be excluded
params.qc_hetfilter  = "both"             // Heterozygosity filter, must be "low", "high" or "both", cf: ./bin/het_check.R

params.pr_window     = 200                // Window size for the pruning, in number of variants
params.pr_step       = 50                 // Window sliding size in number of variants
params.pr_r2         = 0.25               // Pairs of variants with r2 > pr_r2 will be removed

params.saige_trait   = "binary"           // Trait type, must be 'binary' or 'quantitative'
params.gwas_maf      = 0.01               // Variants with a MAF < 'gwas_maf' will be removed for the GWAS analysis
params.rvat_maf      = 0.01               // Variants with a MAF > 'rvat_maf' will be removed for the rare variants analysis


// Checking input values:
if(params.plink_fileset == "")              error("\nERROR in config: 'plink_fileset' is required")

if(params.qc_mind < 0)                      error("\nERROR in config: 'qc_mind' must be >= 0")
if(params.qc_geno < 0)                      error("\nERROR in config: 'qc_geno' must be >= 0")
if(params.qc_hetfilter != "low" && 
   params.qc_hetfilter != "high" &&  
   params.qc_hetfilter != "both")           error("\nERROR in config: 'qc_hetfilter' must be 'low', 'high' or 'both'")

if(params.pr_window < 1)                    error("\nERROR in config: 'pr_window' must be > 0")
if(params.pr_step < 1)                      error("\nERROR in config: 'pr_step' must be > 0")
if(params.pr_r2 < 0 || params.pr_r2 > 1)    error("\nERROR in config: 'pr_r2' must be between 0 and 1")

if(params.rvat_maf < 0)                     error("\nERROR in config: 'rvat_maf' must be >= 0")
if(params.gwas_maf < 0)                     error("\nERROR in config: 'gwas_maf' must be >= 0")

if(params.saige_trait != "binary" &&
   params.saige_trait != "quantitative")    error("\nERROR in config: 'saige_trait' must be 'binary' or 'quantitative'")


// Subworkflows QC & GWAS & RVAT:
workflow QC {
    take:
        plink_ch
        remove_ch

    main:
        BaseQC(plink_ch, remove_ch)
        plink_baseQC_ch = BaseQC.out.plink_baseQC

        Pruning(plink_baseQC_ch)
        HetCoeff(plink_baseQC_ch, Pruning.out.prune_in)
        HetFilter(HetCoeff.out.het)
        CreateOutputGWAS(plink_baseQC_ch, HetFilter.out.valides)
        PlotPCA(CreateOutputGWAS.out.plink_QCed, CreateOutputGWAS.out.eigenvec)

    emit:
        CreateOutputGWAS.out.plink_QCed
}

workflow SAIGE_GWAS {
    take:
        plink_QCed
}

workflow SAIGE_RVAT {

}


workflow {
    plink_ch = Channel.fromFilePairs(params.plink_fileset, size: 3)
    params.qc_remove == "" ? remove_ch = [] : remove_ch = Channel.fromPath(params.qc_remove)
    plink_QCed = QC(plink_ch, remove_ch)
}