#!/usr/bin/env Rscript

# Script to plot a PCA with admixture Pie Charts on top of each point.
if(!require(optparse, quietly = T)) install.packages(optparse) 
if(!require(ggplot2, quietly = T)) install.packages(ggplot2) 

# Function to make the pie plot for each individual:
make_pie_polygons <- function(values, cx, cy, radius = 0.05, id) {
    angles <- c(0, cumsum(values)) * 2 * pi
  
    res <- lapply(seq_along(values), function(i) {
                theta <- seq(angles[i], angles[i+1], length.out = 40)
                x <- cx + radius * cos(theta)
                y <- cy + radius * sin(theta)

                data.frame(pie_id = id,
                           slice  = i,
                           x = c(cx, x),
                           y = c(cy, y)) 
                })
  do.call(rbind, res)
}

# Parsing the option from the Rscript command
option_list <- list(
    make_option(c("-e", "--eigenvec"), 
                type = "character", 
                default = "", 
                help = "Full path to plink .eigenvec file.",
                metavar = "<PATH_TO_EIGENVEC>"),

    make_option(c("-a", "--admixture"),
                type = "character",
                default = "",
                help = "Full path to the admixture table.",
                metavar = "<PATH_TO_ADXM_TABLE>"),

    make_option(c("-o", "--output"),
                type = "character",
                default = "admixture_pie",
                help = "Output file name",
                metavar = "<OUTPUT_BASE>"),

    make_option(c("-v", "--verbose"),
                type = "logical",
                default = FALSE,
                help = "Print information during execution (FALSE by default)",
                metavar = "<TRUE,FALSE>")
);

opt <- parse_args(OptionParser(option_list = option_list))

# Reading the file with the eigenvec data:
if(file.exists(opt$e)) {
    if(opt$v) print(paste0("Reading the .eigenvec file: ", opt$e))
    eigenvec <- read.table(file = opt$e, header = F, sep = " ", stringsAsFactors = F)
    colnames(eigenvec) <- c("FID", "IID", paste0("PC", c(1:(ncol(eigenvec)-2))))
} else { stop(paste0("File '", opt$e, "' does not exist.")) }

# Reading the file with the admixture table:
if(file.exists(opt$a)) {
    if(opt$v) print(paste0("Reading the admixture table: ", opt$a))
    admixture <- read.table(file = opt$a, header = F, sep = " ", stringsAsFactors = F)
    colnames(admixture) <- paste0("Pop", c(1:ncol(admixture)))
} else { stop(paste0("File '", opt$a, "' does not exist.")) }

# Selecting IDD,PC1,PC2 and admixture values:
pca_admx <- cbind(eigenvec, admixture)
pca_admx_wide <- pca_admx[,c(2:4,(ncol(eigenvec)+1):(ncol(eigenvec)+ncol(admixture)))]

pie_polys <- do.call(rbind, lapply(1:nrow(pca_admx_wide), function(i) {
        make_pie_polygons(values = unlist(pca_admx_wide[i, colnames(admixture)]),
                          cx = pca_admx_wide$PC1[i],
                          cy = pca_admx_wide$PC2[i],
                          radius = 0.01,
                          id = pca_admx_wide$IID[i])})
                    )

xrange <- range(pca_admx_wide$PC1)
yrange <- range(pca_admx_wide$PC2)
L <- max(diff(xrange), diff(yrange))
xlim_corrected <- c(mean(xrange) - L/2, mean(xrange) + L/2)
ylim_corrected <- c(mean(yrange) - L/2, mean(yrange) + L/2)

pdf(opt$o, width = 8, height = 8)
    ggplot(pie_polys, aes(x = x, y = y, group = interaction(pie_id, slice), fill = factor(slice))) +
        geom_polygon(color = "black") +
        labs(x = "PC1", y = "PC2", fill = "Population") +
        coord_fixed(ratio = 1) +
        scale_x_continuous(limits = xlim_corrected, expand = c(0, 0)) + # To avoid margins on the x axis
        scale_y_continuous(limits = ylim_corrected, expand = c(0, 0)) + # To avoid margins on the y axis
        theme_bw()
dev.off()
