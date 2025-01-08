# Script to create the phenoFile to input in SAIGE+ for the GWAS and RVAT analyses:

# Takes a .fam file and a .covariates file in the plink format and outputs a file
# corresponding to SAIGE format: id + covars + phenotype.

# Our nextflow expects the column with the ids to be named as "IID" and the column with
# the phenotypes as "PHENOTYPE".

# The covariates column names must be defined in the config file.

if(!require(optparse, quietly = T)) install.packages(optparse) 

option_list <- list(
    make_option(c("-f", "--fam"), 
                type = "character", 
                default = "", 
                help = "Full path to plink .fam file.",
                metavar = "character"),
    
    make_option(c("-c", "--covar"),
                type = "character",
                default = "",
                help = "Full path to plink .covariates file.",
                metavar = "character"),

    make_option(c("-v", "--verbose"),
                type = "logical",
                default = FALSE,
                help = "Print information during execution (off by default)")
);

opt <- parse_args(OptionParser(option_list = option_list))

# Reading the .fam file
if(file.exists(opt$f)) {
    if(opt$v) print(paste0("Reading the .fam file", opt$i))
    fam <- read.table(opt$f, header = F, sep = " ")
    colnames(fam) <- c("FID", "IID", "FATHER", "MOTHER", "SEX", "PHENOTYPE")
} else { stop(paste0("File '", opt$f, "' does not exist.")) }

# Reading the .covar file
if(file.exists(opt$c)) {
    if(opt$v) print(paste0("Reading the .covariates file", opt$i))
    covar <- read.table(opt$c, header = T, sep = "\t")
} else { stop(paste0("File '", opt$c, "' does not exist.")) }

# Adding phenotypes to covar:
covar$PHENOTYPE <- fam$PHENOTYPE[match(x = covar$IID, table = fam$IID)]

write.table(x = covar, file = "saige_phenofile.tsv", quote = F, sep = "\t", row.names = F)
