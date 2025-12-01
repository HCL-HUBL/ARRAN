# Script to make the barplots from Admixture's output:
# Each bar correspond to a sample and the colors corresponds to the % of genome
# coming from each ancestry.

if(!require(optparse, quietly = T)) install.packages(optparse)
if(!require(RColorBrewer, quietly = T)) install.packages(RColorBrewer)

# Parsing the option from the Rscript command
option_list <- list(
    # The .single_variant.tsv file from SAIGE
    make_option(c("-i", "--input"), 
                type = "character", 
                default = "", 
                help = "Path to Admixture .Q file.",
                metavar = "<PATH_TO_Q_FILE>"),

    make_option(c("-o", "--output"),
                type = "character",
                default = "Admixture",
                help = "Name of the output pdf file.",
                metavar = "<OUTPUT_NAME>"),

    make_option(c("-k", "--k"),
                type = "integer",
                default = 2,
                help = "The number of populations that was used during Admixture",
                metavar = "<K>")
);

opt <- parse_args(OptionParser(option_list = option_list))

# Reading the admixture table:
if(file.exists(opt$i)) {
    admixture <- read.table(opt$i, header = F, sep = " ")
} else { stop(paste0("File '", opt$i, "' does not exist.")) }

if(opt$k < 0) stop(paste0("K must be > 0, current value '", opt$k,"'"))
if(opt$k > 7) stop(paste0("K must be < 8, current value '", opt$k,"'"))

admixture$pop <- names(admixture)[max.col(admixture)]

admixture_sorted <- data.frame()
for(current_pop in unique(admixture$pop)) {
    subset <- admixture[admixture$pop == current_pop, -ncol(admixture)]
    subset <- subset[order(subset[,1]),]
    admixture_sorted <- rbind(admixture_sorted, subset)
}

pdf(opt$o, width = 12)
    barplot(height = t(as.matrix(admixture_sorted)),  
            col = RColorBrewer::brewer.pal(name = "Dark2", n = 7)[c(1:opt$k)],
            xlab = "Individuals", ylab = "Ancestry Fraction", border = NA, xaxt = 'n')
dev.off()

