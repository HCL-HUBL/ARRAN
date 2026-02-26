#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.enable.moduleBinaries = true

// Importing Processes from modules:
include { GenotypesPreprocessing }  from './modules/GenotypesPreprocessing'
include { BaseQC }                  from './modules/BaseQC'
include { Pruning }                 from './modules/Pruning'
include { HetCoeff }                from './modules/Het/HetCoeff'
include { HetFilter }               from './modules/Het/HetFilter'
include { HWEFlag }                 from './modules/HWEFlag'
include { CreateEigenvec }          from './modules/PCA/CreateEigenvec'
include { PlotPCA }                 from './modules/PCA/PlotPCA'
include { CreateOutputBaseQC }      from './modules/CreateOutputBaseQC'

include { AdmixtureRun }            from './modules/Admixture/AdmixtureRun'
include { AdmixtureBarplot }        from './modules/Admixture/AdmixtureBarplot'
include { AdmixturePieplot }        from './modules/Admixture/AdmixturePieplot'

// Processes to split the data into autosomes and chrX (and then into GWAS and RVAS)
// Process that will split the SNPs belonging to the autosomes (1-22 + 25(PAR)) vs. chrX (23)
include { SplitAutosomesChrX }      from './modules/Split/SplitAutosomesChrX'
include { CreateOutputGWAS }        from './modules/Split/CreateOutputGWAS'
include { CreateOutputRVAT }        from './modules/Split/CreateOutputRVAT'

include { CreatePhenoFile }         from './modules/SAIGE/CreatePhenoFile'
include { CreateSparseGRM }         from './modules/SAIGE/CreateSparseGRM'
include { SaigeFitNullModel }       from './modules/SAIGE/SaigeFitNullModel'
include { SaigeSingleAssoc }        from './modules/SAIGE/SaigeSingleAssoc'
include { CreateGroupFile }         from './modules/SAIGE/CreateGroupFile'
include { SaigeGeneAssoc }          from './modules/SAIGE/SaigeGeneAssoc'

include { XwasQC }                  from './modules/XWAS/XwasQC'
include { XwasSNVsAssoc }           from './modules/XWAS/XwasSNVsAssoc'

include { QQPlot }                  from './modules/QQPlot'
include { QQPlot as QQPlot_comb }   from './modules/QQPlot'
include { QQPlot as QQPlot_M }      from './modules/QQPlot'
include { QQPlot as QQPlot_F }      from './modules/QQPlot'
include { SummaryStatistics }       from './modules/SummaryStatistics'
include { ManhattanPlot }           from './modules/ManhattanPlot'
include { RegionPlot }              from './modules/RegionPlot'
include { ManhattanPlot as ManhattanPlot_comb }      from './modules/ManhattanPlot'
include { ManhattanPlot as ManhattanPlot_M }         from './modules/ManhattanPlot'
include { ManhattanPlot as ManhattanPlot_F }         from './modules/ManhattanPlot'


// Initialising the options:
// General options:
params.plink_fileset    = ""
params.covar_file       = ""                   // Path to file listing the covariates to be included in the model
params.annot_file       = ""                   // (optional) path to the SNP annotations (format: "SNP annotation", space delimited)
params.genome_build     = "hg19"               // Accepts 'hg19' or 'hg38': used to define the genes list and PAR regions to use
params.trait_type       = "binary"             // Trait type, must be 'binary' or 'quantitative'

// Output options:
params.out_basename     = "arran_"             // Basename of the output files
params.outdir           = "${launchDir}"       // Path to the output folder

// QC options:
params.qc_mind          = 0.05                 // Individuals with >5% missing genotypes will be removed
params.qc_geno          = 0.05                 // Variants with >5% missing genotypes will be removed
params.qc_remove        = ""                   // (optional) File listing IIDs and FIDs of individuals to be excluded
params.qc_hetfilter     = "both"               // Heterozygosity filter, must be "none", "low", "high" or "both", cf: ./bin/het_check.R
params.qc_hwe           = 5e-6                 // Variants with a HWE exact test p-value < 'qc_hwe' will be flagged (added to a file for further inspection) 

// Pruning options:
params.pr_window        = 200                  // Window size for the pruning, in number of variants
params.pr_step          = 50                   // Window sliding size in number of variants
params.pr_r2            = 0.25                 // Pairs of variants with r2 > pr_r2 will be removed

params.makeGRM          = false

// GWAS + XWAS options:
params.run_GWAS         = true                  // Boolean indicating whether to run the GWAS analysis or not
params.run_XWAS         = true                  // Boolean indicating whether to run the XWAS analysis or not
params.gwas_maf         = 0.01                  // Variants with a MAF < 'gwas_maf' will be removed for the GWAS analysis
params.xwas_alpha       = 0.05                  // Significance for the X-specific QC steps (bonferroni correction will be applied to this threshold)
params.xwas_covar       = "PC1,PC2,PC3,PC4,PC5" // Covariates for the chrX analysis with XWAS

// RVAT options:
params.run_RVAT         = true                  // Boolean indicating whether to run the Rare Variants analysis or not
params.saige_annot      = "no_annot"            // How to test the annotations (see annotation_in_groupTest in SAIGE)
params.saige_covar      = "PC1,PC2,PC3,PC4,PC5" // List of all covariates to include in the model, comma separated
params.saige_qcovar     = ""                    // List the covariates which are categorical  
params.saige_regions    = ""                    // (optional) list of regions to analyse (bed format)
params.saige_extension  = 5                     // When assigning SNPs to genes, extends the gene bounds by this many kbp

// Admixture:
params.run_admixture    = false                 // Boolean indicating whether to run Admixture
params.admixture_K      = 2                     // Number of expected populations in the dataset


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
if(params.xwas_alpha < 0)                   error("\nERROR in config: 'xwas_alpha' must be >= 0, current value '${params.xwas_alpha}'")
if(params.admixture_K <= 0)                 error("\nERROR in config: 'admixture_K' must be > 0, current value '${params.admixture_K}'")

// Initialising Channels based on params:
plink_ch        = Channel.fromFilePairs(params.plink_fileset, size: 3)
glist_ch        = Channel.fromPath("${projectDir}/data/glist-${params.genome_build}")
annot_ch        = params.annot_file ? Channel.fromPath(params.annot_file) : []
covar_file_ch   = Channel.fromPath(params.covar_file)
remove_ch       = params.qc_remove ? Channel.fromPath(params.qc_remove) : []
regions_ch      = params.saige_regions ? Channel.fromPath(params.saige_regions) : []


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

        CreateEigenvec(CreateOutputBaseQC.out.plink_QCed)
        PlotPCA(CreateOutputBaseQC.out.plink_QCed, CreateEigenvec.out.eigenvec)

    emit:
        plink_QCed          = CreateOutputBaseQC.out.plink_QCed
        plink_QCed_pruned   = CreateOutputBaseQC.out.plink_QCed_pruned
        eigenvec            = CreateEigenvec.out.eigenvec
}


workflow ADMIXTURE {
    take:
        plink_QCed_pruned
        eigenvec

    main:
        AdmixtureRun(plink_QCed_pruned)
        AdmixtureBarplot(AdmixtureRun.out.admixture_table)
        AdmixturePieplot(AdmixtureRun.out.admixture_table, eigenvec)
}


workflow SAIGE_GWAS {
    take:
        autosomes_QCed
        phenoFile_ch
        regions_ch
        glist

    main:
        CreateOutputGWAS(autosomes_QCed, regions_ch)
        
        SaigeFitNullModel(CreateOutputGWAS.out.plink_GWAS, 
                          phenoFile_ch, 
                          "GWAS")

        SaigeSingleAssoc(CreateOutputGWAS.out.plink_GWAS, 
                         SaigeFitNullModel.out.gmmat, 
                         SaigeFitNullModel.out.vr)

        ManhattanPlot(SaigeSingleAssoc.out.saige_sv, "CHR", "POS", "MarkerID", "p.value")
        QQPlot(SaigeSingleAssoc.out.saige_sv, "p.value")
        SummaryStatistics(SaigeSingleAssoc.out.saige_sv, "SAIGE")
        RegionPlot(CreateOutputGWAS.out.plink_GWAS, SummaryStatistics.out, glist)

    emit:
        SaigeSingleAssoc.out.saige_sv
}


workflow SAIGE_RVAT {
    take:
        autosomes_QCed
        phenoFile_ch
        regions_ch
        glist
        annot
    
    main:
        CreateOutputRVAT(autosomes_QCed, regions_ch)

        SaigeFitNullModel(CreateOutputRVAT.out.plink_RVAT, 
                          phenoFile_ch, 
                          "RVAT")

        CreateGroupFile(CreateOutputRVAT.out.plink_RVAT, glist, annot)

        SaigeGeneAssoc(CreateOutputRVAT.out.plink_RVAT, 
                       SaigeFitNullModel.out.gmmat, 
                       SaigeFitNullModel.out.vr, 
                       CreateGroupFile.out.group_file)

        QQPlot(SaigeGeneAssoc.out.saige_gene, "Pvalue")

    emit:
        SaigeGeneAssoc.out.saige_gene
}


workflow XWAS {
    take:
        chrX_basename
        chrX_bed
        chrX_bim
        chrX_fam
        phenoFile_ch

    main:
        n_var = chrX_bim.countLines()

        // X-specific QC steps only applicable to binary traits:
        if(params.trait_type == "binary") {
            XwasQC(chrX_basename, chrX_bed, chrX_bim, chrX_fam, n_var)
            chrX_ch = XwasQC.out.plink_chrX_QCed
        } else {
            chrX_bed_pairs = chrX_basename.combine(chrX_bed)
            chrX_bim_pairs = chrX_basename.combine(chrX_bim)
            chrX_fam_pairs = chrX_basename.combine(chrX_fam)
            chrX_ch = chrX_bed_pairs.join(chrX_bim_pairs).join(chrX_fam_pairs).map { [ it[0], [ it[1], it[2], it[3] ] ] }
        }

        XwasSNVsAssoc(chrX_ch, phenoFile_ch)

        ManhattanPlot_comb(XwasSNVsAssoc.out.xstrat_assoc, "CHR", "BP", "SNP", "P_comb_Fisher")
        ManhattanPlot_M(XwasSNVsAssoc.out.xstrat_assoc, "CHR", "BP", "SNP", "P_M")
        ManhattanPlot_F(XwasSNVsAssoc.out.xstrat_assoc, "CHR", "BP", "SNP", "P_F")

        QQPlot_comb(XwasSNVsAssoc.out.xstrat_assoc, "P_comb_Fisher")
        QQPlot_M(XwasSNVsAssoc.out.xstrat_assoc, "P_M")
        QQPlot_F(XwasSNVsAssoc.out.xstrat_assoc, "P_F")
}


// Main workflow, calling all the other subworkflow:
workflow {
    // Perform base Quality Control on the genotype data:
    QC(plink_ch, remove_ch)

    if(params.run_admixture) {
        ADMIXTURE(QC.out.plink_QCed_pruned, QC.out.eigenvec)
    }

    // Create the sparse GRM (for SAIGE):
    if(params.makeGRM) {
        CreateSparseGRM(QC.out.plink_QCed_pruned)
    }

    // Split autosomes and chrX (chrX will be subjected to specific QC and association tests):
    SplitAutosomesChrX(QC.out.plink_QCed)
    
    CreatePhenoFile(SplitAutosomesChrX.out.autosomes, // In practice: we only need the .fam file
                    QC.out.eigenvec, 
                    covar_file_ch)

    if(params.run_GWAS) {
        SAIGE_GWAS(SplitAutosomesChrX.out.autosomes,
                   CreatePhenoFile.out.phenoFile, 
                   regions_ch,
                   glist_ch)
    }

    if(params.run_XWAS) {
        XWAS(SplitAutosomesChrX.out.chrX_basename,
             SplitAutosomesChrX.out.chrX_bed,
             SplitAutosomesChrX.out.chrX_bim,
             SplitAutosomesChrX.out.chrX_fam,
             CreatePhenoFile.out.phenoFile)
    }

    if(params.run_RVAT) {
        SAIGE_RVAT(SplitAutosomesChrX.out.autosomes,
                   CreatePhenoFile.out.phenoFile,
                   regions_ch,
                   glist_ch,
                   annot_ch)
    }
}
