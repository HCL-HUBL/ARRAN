# Script to create a manhattan plot from SAIGE single variants association output
# Writes it to a pdf file

if(!require(optparse, quietly = T)) install.packages(optparse) 
if(!require(ggplot2, quietly = T))  install.packages(ggplot2) 
if(!require(qqman, quietly = T))    install.packages(qqman) 

# Parsing the option from the Rscript command
option_list <- list(
    # The .single_variant.tsv file from SAIGE
    make_option(c("-i", "--input"), 
                type = "character", 
                default = "", 
                help = "Full path to SAIGE .single_variant.tsv file.",
                metavar = "character"),

    make_option(c("-o", "--output"),
                type = "character",
                default = "Manhattan",
                help = "Name of the output pdf file.",
                metavar = "character"),

    make_option(c("-v", "--verbose"),
                type = "logical",
                default = FALSE,
                help = "Print information during execution (FALSE by default)",
                metavar = "BOOLEAN")
);

opt <- parse_args(OptionParser(option_list = option_list))

# Reading the .het file
if(file.exists(opt$i)) {
    if(opt$v) print(paste0("Reading the .single_variant.tsv file", opt$i))
    saige_sv <- read.table(opt$i, header = T, sep = "\t")
} else { stop(paste0("File '", opt$i, "' does not exist.")) }

pdf(opt$o, width = 12)
    qqman::manhattan(x = saige_sv,
                     chr = "CHR", 
                     bp = "POS", 
                     snp = "MarkerID", 
                     p = "p.value", 
                     col = c("#4393c3", "#009988"), 
                     main = "Manhattan Plot")
dev.off()
