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
include { RunAdmixture }            from './modules/QC.nf'
include { PlotAdmixture }           from './modules/QC.nf'

include { Split_Autosomes_ChrX }    from './modules/Split.nf'
include { CreateOutputGWAS }        from './modules/Split.nf'
include { CreateOutputRVAT }        from './modules/Split.nf'

include { CreatePhenoFile }         from './modules/SAIGE.nf'
include { CreateSparseGRM }         from './modules/SAIGE.nf'
include { SaigeFitNullModel }       from './modules/SAIGE.nf'
include { SaigeSingleAssoc }        from './modules/SAIGE.nf'
include { CreateGroupFile }         from './modules/SAIGE.nf'
include { SaigeGeneAssoc }          from './modules/SAIGE.nf'

include { ChrX_specific_QC }        from './modules/XWAS.nf'
include { ChrX_SNVs_Assoc }         from './modules/XWAS.nf'

include { ManhattanPlot }           from './modules/Downstream.nf'
include { QQPlot }                  from './modules/Downstream.nf'


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
params.rvat_maf         = 0.01                 // Variants with a MAF > 'rvat_maf' will be removed for the RVAT analysis

// Admixture:
params.admixture_K      = 2                    // Number of expected populations in the dataset

// SAIGE options:
params.saige_covar      = "AGE,SEX,PC1,PC2"    // List of all covariates to include in the model, comma separated
params.saige_qcovar     = "SEX"                // List the covariates which are categorical  
params.saige_regions    = ""                   // (optional) list of regions to analyse (bed format)
params.saige_extension  = 5                    // When assigning SNPs to genes, extends the gene bounds by this many kbp

// XWAS options:
params.alpha            = 0.05                 // Significance for the X-specific QC steps (bonferroni correction will be applied to this threshold)

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
if(params.rvat_maf < 0)                     error("\nERROR in config: 'rvat_maf' must be <= 0, current value '${params.rvat_maf}'")

if(params.admixture_K <= 0)                 error("\nERROR in config: 'admixture_K' must be > 0, current value '${params.admixture_K}'")

if(params.xwas_alpha < 0)                   error("\nERROR in config: 'xwas_alpha' must be >= 0, current value '${params.xwas_alpha}'")

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

    main:
        GenotypesPreprocessing(plink_ch)
        BaseQC(GenotypesPreprocessing.out.plink_preprocessed, remove_ch)
        Pruning(BaseQC.out.plink_baseQC)
        HetCoeff(BaseQC.out.plink_baseQC, Pruning.out.prune_in)
        HetFilter(HetCoeff.out.het)
        
        CreateOutputBaseQC(BaseQC.out.plink_baseQC, Pruning.out.prune_in, HetFilter.out.valides) // BaseQC output
        
        HWEFlag(CreateOutputBaseQC.out.plink_QCed, params.qc_hwe)

        CreateEigenvec(CreateOutputBaseQC.out.plink_QCed_pruned)
        PlotPCA(CreateOutputBaseQC.out.plink_QCed_pruned, CreateEigenvec.out.eigenvec)

        // RunAdmixture(CreateOutputBaseQC.out.plink_QCed_pruned)
        // PlotAdmixture(RunAdmixture.out.admixture_table)

    emit:
        plink_QCed          = CreateOutputBaseQC.out.plink_QCed
        plink_QCed_pruned   = CreateOutputBaseQC.out.plink_QCed_pruned
        eigenvec            = CreateEigenvec.out.eigenvec
}


// Subworklow for the common variant Single Association analysis with SAIGE+
workflow SAIGE_GWAS {
    take:
        autosomes_QCed
        sparse_GRM
        sparse_ids
        phenoFile_ch
        regions_ch

    main:
        CreateOutputGWAS(autosomes_QCed, regions_ch)
        
        SaigeFitNullModel(CreateOutputGWAS.out.plink_GWAS, 
                          sparse_GRM, sparse_ids, 
                          phenoFile_ch, 
                          "GWAS")

        SaigeSingleAssoc(CreateOutputGWAS.out.plink_GWAS, 
                         sparse_GRM, sparse_ids, 
                         SaigeFitNullModel.out.gmmat, 
                         SaigeFitNullModel.out.vr)

        ManhattanPlot(SaigeSingleAssoc.out.saige_sv, "CHR", "POS", "MarkerID", "p.value")
        QQPlot(SaigeSingleAssoc.out.saige_sv, "p.value")

    emit:
        SaigeSingleAssoc.out.saige_sv
}


// Subworklow for the Rare Variant Association Test with SAIGE+
workflow SAIGE_RVAT {
    take:
        autosomes_QCed
        sparse_GRM
        sparse_ids
        phenoFile_ch
        regions_ch
        glist
    
    main:
        CreateOutputRVAT(autosomes_QCed, regions_ch)

        SaigeFitNullModel(CreateOutputRVAT.out.plink_RVAT, 
                          sparse_GRM, sparse_ids, 
                          phenoFile_ch, 
                          "RVAT")

        CreateGroupFile(CreateOutputRVAT.out.plink_RVAT, glist)

        SaigeGeneAssoc(CreateOutputRVAT.out.plink_RVAT, 
                       sparse_GRM, sparse_ids, 
                       SaigeFitNullModel.out.gmmat, 
                       SaigeFitNullModel.out.vr, 
                       CreateGroupFile.out.group_file)

        QQPlot(SaigeGeneAssoc.out.saige_gene, "Pvalue")

    emit:
        SaigeGeneAssoc.out.saige_gene
}


workflow XWAS {
    take:
        plink_chrX_basename
        plink_chrX_bed
        plink_chrX_bim
        plink_chrX_fam
        phenoFile_ch

    main:
        n_var = plink_chrX_bim.countLines()

        if(params.trait_type == "binary") {
            ChrX_specific_QC(plink_chrX_basename, plink_chrX_bed, plink_chrX_bim, plink_chrX_fam, n_var)
            chrX_ch = ChrX_specific_QC.out.plink_chrX_QCed
        } else {
            chrX_ch = plink_chrX
        }

        ChrX_SNVs_Assoc(chrX_ch, phenoFile_ch)

        ManhattanPlot(ChrX_SNVs_Assoc.out.xstrat_assoc, "CHR", "BP", "SNP", "P_F")
}


// Main workflow, calling all the other subworkflow:
workflow {

    // Perform base Quality Control on the genotype data:
    QC(plink_ch, remove_ch)

    // Create the sparse GRM (for SAIGE):
    CreateSparseGRM(QC.out.plink_QCed_pruned)

    // Split autosomes and chrX (chrX will be subjected to specific QC and association tests):
    Split_Autosomes_ChrX(QC.out.plink_QCed)
    
    // Run SAIGE+ (GWAS and RVAT) on the autosomes:
    CreatePhenoFile(Split_Autosomes_ChrX.out.autosomes,  // In practice: we only need the .fam file not the complete plink fileset
                    QC.out.eigenvec, 
                    covar_file_ch)

    SAIGE_GWAS(Split_Autosomes_ChrX.out.autosomes,
               CreateSparseGRM.out.sparseGRM,
               CreateSparseGRM.out.sampleIDs,
               CreatePhenoFile.out.phenoFile, 
               regions_ch)

    SAIGE_RVAT(Split_Autosomes_ChrX.out.autosomes,
               CreateSparseGRM.out.sparseGRM,
               CreateSparseGRM.out.sampleIDs,
               CreatePhenoFile.out.phenoFile,
               regions_ch,
               glist_ch)

    // Run XWAS (chrX-specific QC, GWAS and RVAT) on chrX:
    XWAS(Split_Autosomes_ChrX.out.chrX_basename,
         Split_Autosomes_ChrX.out.chrX_bed,
         Split_Autosomes_ChrX.out.chrX_bim,
         Split_Autosomes_ChrX.out.chrX_fam,
         CreatePhenoFile.out.phenoFile)
}
