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
    if(opt$v) print(paste0("Reading the .eigenvec file: ", opt$i))
    eigenvec <- read.table(file = opt$i, header = F, sep = " ", stringsAsFactors = F)
    colnames(eigenvec) <- c("FID", "IID", paste0("PC", c(1:(ncol(eigenvec)-2))))
} else { stop(paste0("File '", opt$i, "' does not exist.")) }

# If a .fam file is provided, we will color the PCA with the phenotype information:
if(file.exists(opt$f)) {
    if(opt$v) print(paste0("Reading the .fam file: ", opt$f))
    fam <- read.table(file = opt$f, header = F, sep = " ", stringsAsFactors = F)
    colnames(fam) <- c("FID", "IID", "FATHER", "MOTHER", "SEX", "PHENO")
    # Reordering the fam file according to the eigenvec, to extract the phenotypes:
    eigenvec$phenotype <- as.factor(fam$PHENO[match(x = eigenvec$IID, table = fam$IID)])
    eigenvec$sex <- as.factor(fam$SEX[match(x = eigenvec$IID, table = fam$IID)])
} else { 
    eigenvec$phenotype <- "no_phenotype_info"
    eigenvec$sex <- "no_sex_info"
}

gg_pca <- ggplot(eigenvec, aes(x = PC1, y = PC2, col = phenotype)) + 
            geom_point() + theme_bw()

gg_pca_sex <- ggplot(eigenvec, aes(x = PC1, y = PC2, col = sex)) + 
                geom_point() + theme_bw()

pdf(paste0(opt$o,".pdf"))
    print(gg_pca)
    print(gg_pca_sex)
dev.off()