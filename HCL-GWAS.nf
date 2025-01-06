#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BaseQC }          from './modules/QC.nf'
include { Pruning }         from './modules/QC.nf'
include { HetCoeff }        from './modules/QC.nf'
include { HetFilter }       from './modules/QC.nf'

// Initialising the options with default values:
params.plink_fileset = ""
params.outdir = "${launchDir}"
params.genome_build = "hg19"

// QC options:
params.qc_mind = 0.05                        // Individuals with >5% missing genotypes will be removed
params.qc_geno = 0.05                        // Variants with >5% missing genotypes will be removed
params.qc_maf  = 0.01                        // Variants with a MAF <1% will be removed 
params.qc_remove = ""                        // (optional) File listing IIDs and FIDs of individuals to be excluded

params.qc_hetfilter = "both"                   // Heterozygosity filter, must be "low", "high" or "both", cf: ./bin/het_check.R

// Pruning options:
params.pr_window = 200                       // Window size for the pruning, in number of variants
params.pr_step = 50                          // Window sliding size in number of variants
params.pr_r2 = 0.25                          // Pairs of variants with r2 > pr_r2 will be removedparams.

// Check input values:
if(params.qc_mind < 0)    error("\nERROR in config: qc_mind must be >= 0")
if(params.qc_geno < 0)    error("\nERROR in config: qc_geno must be >= 0")
if(params.qc_maf < 0)     error("\nERROR in config: qc_mad must be >= 0")

if(params.qc_hetfilter != "low" && 
    params.qc_hetfilter != "high" && 
     params.qc_hetfilter != "both")      error("\nERROR in config: 'qc_filter' must be 'low', 'high' or 'both'")


workflow {
    println("Working directory: '"+params.outdir+"'")

    geno_ch = Channel.fromFilePairs(params.plink_fileset, size: 3)

    if(params.qc_remove  == "") { remove_ch = [] } else { remove_ch = Channel.fromPath(params.qc_remove) }

    BaseQC(geno_ch, remove_ch)
    plink_QCed_ch = BaseQC.out.plink_QCed
    eigenvec_ch   = BaseQC.out.eigenvec
    
    Pruning(plink_QCed_ch)
    prune_in_ch = Pruning.out.prune_in
    
    HetCoeff(plink_QCed_ch, prune_in_ch)

    HetFilter(HetCoeff.out.het)
    valides_ch = HetFilter.out.valides
}