# Script to plot a PCA based on a .eigenvec file produced with plink (option --pca).

# This script plots the PCA in pdf and png files.

if(!require(optparse, quietly = T)) install.packages(optparse) 
if(!require(ggplot2, quietly = T)) install.packages(ggplot2) 

# Parsing the option from the Rscript command
option_list <- list(
    make_option(c("-i", "--input"), 
                type = "character", 
                default = "", 
                help = "Full path to plink .eigenvec file.",
                metavar = "character"),

    make_option(c("-f", "--fam"),
                type = "character",
                default = "",
                help = "Optional: the path to plink .fam file to color the PCA.",
                metavar = "character"),
    
    make_option(c("-o", "--output"),
                type = "character",
                default = "PCA",
                help = "Output basename, the PCA will be plotted to <output>.pdf and <output>.png",
                metavar = "character"),

    make_option(c("-v", "--verbose"),
                type = "logical",
                default = FALSE,
                help = "Print information during execution (off by default)")
);

opt <- parse_args(OptionParser(option_list = option_list))

# Reading the file with the eigenvec data:
if(file.exists(opt$i)) {
    if(opt$v) print(paste0("Reading .eigenvec file: ", opt$i))
    eigenvec <- read.table(file = opt$i, header = F, sep = " ", stringsAsFactors = F)
    colnames(eigenvec_mx) <- c("FID", "IID", paste0("PC", c(1:ncol(eigenvec_mx))))
} else { stop(paste0("File '", opt$i, "' does not exist.")) }

# If a .fam file is provided, we will color the PCA with the phenotype information:
if(file.exists(opt$f)) {
    if(opt$v) print(paste0("Reading .fam file: ", opt$f))
    fam <- read.table(file = opt$f, header = T, sep = " ", stringsAsFactors = F)

    eigenvec$phenotype <-  
} else { eigenvec$phenotype <- "no_phenotype" }





gg_pca <- ggplot(eigenvec_mx, aes(x = PC1, y = PC2, col = design)) + 
            geom_point() + theme_bw()
