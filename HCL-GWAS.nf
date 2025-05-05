#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Importing Processes from modules:
include { GenotypesPreprocessing }  from './modules/QC.nf'
include { BaseQC }                  from './modules/QC.nf'
include { PlotPCA }                 from './modules/QC.nf'
include { Pruning }                 from './modules/QC.nf'
include { HetCoeff }                from './modules/QC.nf'
include { HetFilter }               from './modules/QC.nf'
include { HWEFlag }                 from './modules/QC.nf'
include { CreateEigenvec }          from './modules/QC.nf'
include { CreateOutputBaseQC }      from './modules/QC.nf'
include { CreateOutputGWAS }        from './modules/QC.nf'
include { CreateOutputRVAT }        from './modules/QC.nf'
include { RunAdmixture }            from './modules/QC.nf'
include { PlotAdmixture }           from './modules/QC.nf'

include { CreatePhenoFile }         from './modules/Association.nf'
include { CreateSparseGRM }         from './modules/Association.nf'
include { SaigeFitNullModel }       from './modules/Association.nf'
include { SaigeSingleAssoc }        from './modules/Association.nf'
include { CreateGroupFile }         from './modules/Association.nf'
include { SaigeGeneAssoc }          from './modules/Association.nf'
include { ManhattanPlot }           from './modules/Association.nf'
include { QQPlot }                  from './modules/Association.nf'


// Initialising the options with default values:
// General options:
params.plink_fileset    = ""                   // The path to the plink fileset (/path/to/example.{bim,bed,fam})
params.covar_file       = ""
params.outdir           = "${launchDir}"       // The output directory where the outputs will be stored
params.out_basename     = "hcl_gwas"           // Basename of the output files
params.genome_build     = "hg19"               // Accepts 'hg19' or 'hg38', used to define PAR regions
params.trait_type       = "binary"             // Trait type, must be 'binary' or 'quantitative'

 // QC options:
params.qc_mind          = 0.05                 // Individuals with >5% missing genotypes will be removed
params.qc_geno          = 0.05                 // Variants with >5% missing genotypes will be removed
params.qc_remove        = ""                   // (optional) File listing IIDs and FIDs of individuals to be excluded
params.qc_hetfilter     = "both"               // Heterozygosity filter, must be "none", "low", "high" or "both", cf: ./bin/het_check.R
params.qc_hwe           = 5e-6

// Pruning options:
params.pr_window        = 200                  // Window size for the pruning, in number of variants
params.pr_step          = 50                   // Window sliding size in number of variants
params.pr_r2            = 0.25                 // Pairs of variants with r2 > pr_r2 will be removed

// MAF options:
params.gwas_maf         = 0.01                 // Variants with a MAF < 'gwas_maf' will be removed for the GWAS analysis

// Admixture:
params.admixture_K      = 2                    // Number of expected populations in the dataset

// SAIGE options
params.saige_covar      = "AGE,SEX,PC1,PC2"    // List of all covariates to include in the model, comma separated
params.saige_qcovar     = "SEX"                // List the covariates which are categorical  
params.saige_regions    = ""                   // (optional) list of regions to analyse (bed format)
params.saige_extension  = 5                    // When assigning SNPs to genes, extends the gene bounds by this many kbp


// Checking input values:
if(params.plink_fileset == "")              error("\nERROR in config: 'plink_fileset' is required")
if(params.covar_file == "")                 error("\nERROR in config: 'covar_file' is required")

if(params.genome_build != "hg19" &&
   params.genome_build != "hg38")           error("\nERROR in config: 'genome_build' must be 'hg19' or 'hg38', current value '${params.genome_build}'")

if(params.trait_type != "binary" &&
   params.trait_type != "quantitative")     error("\nERROR in config: 'trait_type' must be 'binary' or 'quantitative', current value '${params.trait_type}'")

if(params.qc_mind < 0)                      error("\nERROR in config: 'qc_mind' must be >= 0, current value '${params.qc_mind}'")
if(params.qc_geno < 0)                      error("\nERROR in config: 'qc_geno' must be >= 0, current value '${params.qc_geno}'")

if(params.qc_hetfilter != "none" &&
   params.qc_hetfilter != "low"  && 
   params.qc_hetfilter != "high" &&  
   params.qc_hetfilter != "both")           error("\nERROR in config: 'qc_hetfilter' must be 'none', 'low', 'high' or 'both', , current value '${params.qc_hetfilter}'")

if(params.qc_hwe > 1)                       error("\nERROR in config: 'qc_hwe' must be a <= 1, current value 'qc_hwe'")

if(params.pr_window < 1)                    error("\nERROR in config: 'pr_window' must be > 0, , current value '${params.pr_window}'")
if(params.pr_step < 1)                      error("\nERROR in config: 'pr_step' must be > 0, current value '${params.pr_step}'")
if(params.pr_r2 < 0 || params.pr_r2 > 1)    error("\nERROR in config: 'pr_r2' must be between 0 and 1, current value '${params.pr_r2}'")

if(params.gwas_maf < 0)                     error("\nERROR in config: 'gwas_maf' must be >= 0, current value '${params.gwas_maf}'")
if(params.admixture_K <= 0)                 error("\nERROR in config: 'admixture_K' must be > 0, current value '${params.admixture_K}'")


// Initialising Channels based on params:
plink_ch        = Channel.fromFilePairs(params.plink_fileset, size: 3)
glist_ch        = Channel.fromPath("${projectDir}/data/glist-${params.genome_build}")
covar_file_ch   = Channel.fromPath(params.covar_file)
remove_ch       = params.qc_remove ? Channel.fromPath(params.qc_remove) : []
regions_ch      = params.saige_regions ? Channel.fromPath(params.saige_regions) : []


// Subworklow for the Quality Control on the genotype data
workflow QC {
    take:
        plink_ch
        remove_ch
        regions_ch

    main:
        GenotypesPreprocessing(plink_ch)
        BaseQC(GenotypesPreprocessing.out.plink_preprocessed, remove_ch)
        Pruning(BaseQC.out.plink_baseQC)
        HetCoeff(BaseQC.out.plink_baseQC, Pruning.out.prune_in)
        HetFilter(HetCoeff.out.het)
        CreateOutputBaseQC(BaseQC.out.plink_baseQC, Pruning.out.prune_in, HetFilter.out.valides) // Output for the RVAT (high MAF variants will be filtered during the association test)
        
        HWEFlag(CreateOutputBaseQC.out.plink_QCed, params.qc_hwe)
        
        CreateSparseGRM(CreateOutputBaseQC.out.plink_QCed_pruned)

        CreateEigenvec(CreateOutputBaseQC.out.plink_QCed)
        PlotPCA(CreateOutputBaseQC.out.plink_QCed, CreateEigenvec.out.eigenvec) // PCA on GWAS output (need to remove very low maf variants to avoid errors)

        SplitX(CreateOutputBaseQC.out.plink_QCed)

        CreateOutputGWAS(CreateOutputBaseQC.out.plink_QCed, regions_ch)
        CreateOutputRVAT(CreateOutputBaseQC.out.plink_QCed, regions_ch)

        RunAdmixture(CreateOutputBaseQC.out.plink_QCed_pruned)
        PlotAdmixture(RunAdmixture.out.admixture_table)

    emit:
        plink_QCed_gwas = CreateOutputGWAS.out.plink_GWAS
        plink_QCed_rvat = CreateOutputRVAT.out.plink_RVAT
        eigenvec        = CreateEigenvec.out.eigenvec
}


// Subworklow for the common variant Single Association analysis with SAIGE+
workflow SAIGE_GWAS {
    take:
        plink_QCed
        pheno_file
    
    main:
        SaigeFitNullModel(plink_QCed, pheno_file, "GWAS")
        SaigeSingleAssoc(plink_QCed, SaigeFitNullModel.out.gmmat, SaigeFitNullModel.out.vr)
        ManhattanPlot(SaigeSingleAssoc.out.saige_sv)
        QQPlot(SaigeSingleAssoc.out.saige_sv, "p.value")

    emit:
        SaigeSingleAssoc.out.saige_sv
}


// Subworklow for the Rare Variant Association Test with SAIGE+
workflow SAIGE_RVAT {
    take:
        plink_QCed
        pheno_file
        glist
    
    main:
        SaigeFitNullModel(plink_QCed, pheno_file, "RVAT")
        CreateGroupFile(plink_QCed, glist)
        SaigeGeneAssoc(plink_QCed, SaigeFitNullModel.out.gmmat, SaigeFitNullModel.out.vr, CreateGroupFile.out.group_file)
        QQPlot(SaigeGeneAssoc.out.saige_gene, "Pvalue")

    emit:
        SaigeGeneAssoc.out.saige_gene
}


// Main workflow, calling all the other subworkflow:
workflow {
    // Perform base Quality Control on the genotype data:
    QC(plink_ch, remove_ch, regions_ch)

    // Extract the QCed genotypes, one for the GWAS (with common variants), one for the RVAT (with rare variants):
    plink_QCed_gwas = QC.out.plink_QCed_gwas
    plink_QCed_rvat = QC.out.plink_QCed_rvat
    eigenvec        = QC.out.eigenvec

    // Create the phenotype file and run SAIGE+ (GWAS and RVAT):
    pheno_file_ch = CreatePhenoFile(plink_QCed_gwas, eigenvec, covar_file_ch) // only needs the .fam file, so it should not matter if we use the gwas or rvat plinks
    
    SAIGE_GWAS(plink_QCed_gwas, pheno_file_ch)
    SAIGE_RVAT(plink_QCed_rvat, pheno_file_ch, glist_ch)
}
