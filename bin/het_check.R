# Script to filter individuals with too high or too low heterozygosity (F coeff).
# Samples are filtered if their F coeff is  or  3 SDs from the population mean.

# The heterozygosity coefficient is obtained with plink --het

# This script outputs two files, one with the valid individuals, that passed the
# thresholds, the other with the invalid individuals.The name of the output is 
# based on the input (adding '.valides' and '.nonvalides')

if(!require(optparse, quietly = T)) install.packages(optparse) 
if(!require(ggplot2, quietly = T)) install.packages(ggplot2) 

# Parsing the option from the Rscript command
option_list <- list(
    # The .het file from plink
    make_option(c("-i", "--input"), 
                type = "character", 
                default = "", 
                help = "Full path to plink .het file.",
                metavar = "<PATH_TO_HET>"),
    
    # Filtering on F coeff (! low het. = high F coeff and inversely)
    make_option(c("-f", "--filter"),
                type = "character",
                default = "both",
                help = "Heterozygosity filter to apply, 'none', 'low', 'high' or 'both'",
                metavar = "<none,low,high,both>"),

    make_option(c("-v", "--verbose"),
                type = "logical",
                default = FALSE,
                help = "Print information during execution (FALSE by default)",
                metavar = "<TRUE,FALSE>")
);

opt <- parse_args(OptionParser(option_list = option_list))

# Reading the .het file
if(file.exists(opt$i)) {
    if(opt$v) print(paste0("Reading the .het file", opt$i))
    het <- read.table(opt$i, header = T)
} else { stop(paste0("File '", opt$i, "' does not exist.")) }

# Computing the F coeff mean and sd
het.moy <- mean(het$F)
het.sd <- sd(het$F)
if(opt$v) print(paste0("F coeff mean:  ", het.moy))
if(opt$v) print(paste0("F coeff SD:    ", het.sd, "\n"))

# ! Heterozygosity and F coeff are inversed
# F coeff = 1 - (obs het  expected het)
het$filter <- "valid"
het$filter[het$F > mean(het$F) + 3*sd(het$F)] <- "low_het"
het$filter[het$F < mean(het$F) - 3*sd(het$F)] <- "high_het"

# Get the valid and nonvalid samples (3 SDs from the group mean)
if(opt$f == "none") {
    valides     <- het
    non.valides <- colnames(het)
} else if(opt$f == "high") {
    valides     <- subset(het, filter == "valid" | filter == "low_het")
    non.valides <- subset(het, filter == "high_het")
} else if(opt$f == "low") {
    valides     <- subset(het, filter == "valid" | filter == "high_het")
    non.valides <- subset(het, filter == "low_het")
} else if(opt$f == "both") {
    valides     <- subset(het, filter == "valid")
    non.valides <- subset(het, filter == "low_het" | filter == "high_het")
} else {
    stop(paste0("Option -f / --filter should either be 'none', 'low', 'high' of 'both', current value '", opt$f, "'"))
}

pdf(paste0(opt$i, '.pdf'))
    ggplot(het, aes(x = F, fill = filter)) +
        geom_histogram(col = "black", bins = 50) + 
        scale_fill_manual(values = c("#4393C3", "#B2182B", "#009988")) + 
        theme_bw() + ggtitle("GWAS QC - Heterozygosity")
dev.off()

valides_filename <- paste0(opt$i, '.valides')
if(opt$v) print(paste0("Writing valid samples to ", valides_filename))
write.table(x = valides[,c(1,2)], file = valides_filename, quote = F, row.names = F)

nonvalides_filename <- paste0(opt$i, '.nonvalides')
if(opt$v) print(paste0("Writing nonvalid samples to ", nonvalides_filename))
write.table(x = non.valides[,c(1,2)], file = nonvalides_filename, quote = F, row.names = F)
