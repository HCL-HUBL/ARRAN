# Script to create a QQ plot from a tsv file containing p-valuee
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
                help = "Full path to the .tsv file.",
                metavar = "<PATH_TO_SAIGE>"),

    make_option(c("-p", "--pcol"),
                type = "character",
                default = "p.value",
                help = "Name of the p value column in the input file",
                metavar = "<NAME P column>"),

    make_option(c("-o", "--output"),
                type = "character",
                default = "Manhattan",
                help = "Name of the output pdf file.",
                metavar = "<OUTPUT_NAME>"),

    make_option(c("-v", "--verbose"),
                type = "logical",
                default = FALSE,
                help = "Print information during execution (FALSE by default)",
                metavar = "<TRUE,FALSE>")
);

opt <- parse_args(OptionParser(option_list = option_list))

# Reading the .het file
if(file.exists(opt$i)) {
    if(opt$v) print(paste0("Reading the .tsv file", opt$i))
    saige <- read.table(opt$i, header = T)
} else { stop(paste0("File '", opt$i, "' does not exist.")) }

if(opt$p %in% colnames(saige)) {
    pcolname <- opt$p
} else { stop(paste0("Column '", pcolname, "' not found in the .tsv file."))}

png(opt$o, height = 14, width = 14, res = 300, unit = "cm")
    qqman::qq(saige[,pcolname])
dev.off()
