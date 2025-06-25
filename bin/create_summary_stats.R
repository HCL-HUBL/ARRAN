# Script to create Summary Statistics from SAIGE / XWAS output

# The format is described in: https://github.com/EBISPOT/gwas-summary-statistics-standard
# (from: https://www.ebi.ac.uk/gwas/docs/summary-statistics-format)

if(!require(optparse, quietly = T)) install.packages(optparse) 

option_list <- list(
    make_option(c("-i", "--input"),
                type = "character",
                default = "",
                help = "Full path to the association file",
                metavar = "<PATH_TO_ASSOC>"),

    make_option(c("-s", "--software"),
                type = "character",
                default = "SAIGE",
                help = "Tool used to generated the assoc file",
                metavar = "<SAIGE/XWAS>"),

    make_option(c("-t", "--trait"),
                type = "character",
                default = "binary",
                help = "Trait type, must be 'binary' or 'quantitative'",
                metavar = "<binary/quantitative>"),
    
    make_option(c("-o", "--out"),
                type = "character",
                default = "",
                help = "Name of the output file.",
                metavar = "<OUTFILE>")
);

opt <- parse_args(OptionParser(option_list = option_list))

# Reading the assoc file
if(file.exists(opt$i)) {
    assoc <- read.table(opt$i, header = T, sep = "\t")
} else { stop(paste0("File '", opt$i, "' does not exist.")) }

es_type = ""
if(opt$t == "binary") {
    es_type = "odds_ratio"
} else if(opt$t == "quantitative") {
    es_type = "binary"
} else {
    stop(paste0("Trait type must be 'binary' or 'quantitative', current value: '", opt$t, "'"))
}

# Changing the column names to stay consistent regardless of the software:
if(opt$s == "SAIGE") {
    colnames(assoc) <- c("chromosome", "base_pair_location", "rsid", "other_allele", "effect_allele",
                        "AC_Allele2", "effect_allele_frequency", "MissingRate", es_type, "standard_error", 
                        "Tstat", "var", "p_value", "p.value.NA", "Is.SPA", "AF_case", "AF_ctrl", 
                        "N_case", "N_ctrl")
} else if(opt$s == "XWAS") {
    # colnames(assoc) <- c("chromosome", "rsid", "base_pair_location", "other_allele", "affect_allele",
    #                      "TEST", "OR_M", "P_M", es_type, "P_F", "Fisher_Chi_Squared", "p_value")
    print("XWAS not yet supported due to missing columns in the output")
} else {
    stop(paste0("'", opt$s, "' is not a valid tool name (valid names: 'SAIGE', 'XWAS')"))
}

summary_stats <- assoc[,c("chromosome", "base_pair_location", "effect_allele",
                          "other_allele", es_type, "standard_error",
                          "effect_allele_frequency", "p_value", "rsid")]

write.table(x = summary_stats, file = opt$o, quote = F, sep = "\t", row.names = F)
