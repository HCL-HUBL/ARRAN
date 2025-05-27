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
                help = "Full path to the association result .tsv file.",
                metavar = "<PATH_TO_SAIGE>"),

    make_option(c("-c", "--chr_col"),
                type = "character",
                default = "CHR",
                help = "Name of the chromosome column",
                metavar = "<CHR_col>"),
    
    make_option(c("-b", "--bp_col"),
                type = "character",
                default = "BP",
                help = "Name of the position column",
                metavar = "<BP_col>"),

    make_option(c("-m", "--marker_col"),
                type = "character",
                default = "MarkerID",
                help = "Name of the SNP id column",
                metavar = "<MarkerID_col>"),
    
    make_option(c("-p", "--pvalue_col"),
                type = "character",
                default = "p.value",
                help = "Name of the P_value column",
                metavar = "<PVAL_col>"),

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
    if(opt$v) print(paste0("Reading the .single_variant.tsv file", opt$i))
    saige_sv <- read.table(opt$i, header = T)
} else { stop(paste0("File '", opt$i, "' does not exist.")) }

pval_col <- which(colnames(saige_sv) == opt$p)

pdf(opt$o, width = 14)
    qqman::manhattan(x = saige_sv[!is.na(saige_sv[,pval_col]),],
                     chr = opt$c,
                     bp = opt$b, 
                     snp = opt$m,
                     p = opt$p,
                     suggestiveline = -log10( 1 / nrow(saige_sv)),
                     genomewideline = -log10( 0.05 / nrow(saige_sv)),
                     col = c("#77AADD", "#EE8866"),
                     main = "Manhattan Plot",
                     cex = 0.6, cex.axis = 0.9)
dev.off()
