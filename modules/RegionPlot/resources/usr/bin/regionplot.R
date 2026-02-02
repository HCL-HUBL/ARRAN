#!/usr/bin/env Rscript

# Script for the regional plot, based on the summary statistic file generated with ARRAN
# But should be compatible with any summary statistic file.

if(!require(optparse, quietly = T)) install.packages(optparse)
if(!require(ggplot2, quietly = T)) install.packages(ggplot2)

# Parsing the option from the Rscript command
option_list <- list(
    make_option(c("-i", "--input"), 
                type = "character", 
                default = "", 
                help = "Full path to the summary statistics file.",
                metavar = "<PATH_TO_STATS>"),

    make_option(c("-s", "--snp"),
                type = "character",
                default = "",
                help = "Identifier of the top SNP. The regional plot will be centered around that SNP.",
                metavar = "<RSID>"),
    
    make_option(c("-w", "--winsize"),
                type = "numeric",
                default = 100000,
                help = "Size of the region to plot",
                metavar = "<SIZE>"),
    
    make_option(c("-l", "--ld"),
                type = "character",
                default = "",
                help = "Full path to the LD file.",
                metavar = "<LD_FILE>"),

    make_option(c("-g", "--genelist"),
                type = "character",
                default = "",
                help = "Full path to the file listing gene positions.",
                metavar = "<GLIST-HG38>"),

    make_option(c("-o", "--output"),
                type = "character",
                default = "regionplot.pdf",
                help = "Name of the output file.",
                metavar = "<OUTNAME>"),

    make_option(c("-v", "--verbose"),
                type = "logical",
                default = FALSE,
                help = "Print information during execution (FALSE by default)",
                metavar = "<TRUE,FALSE>")
);

opt <- parse_args(OptionParser(option_list = option_list))

# Reading the file with the summary statistics:
if(file.exists(opt$i)) {
    if(opt$v) print(paste0("Reading the summary statistics file: ", opt$i))
    summ_stat <- read.table(file = opt$i, header = T, sep = "\t", stringsAsFactors = F)
} else { stop(paste0("File '", opt$i, "' does not exist.")) }

# Reading the file with the LD info:
if(file.exists(opt$l)) {
    if(opt$v) print(paste0("Reading the LD file: ", opt$l))
    LD_info <- read.table(file = opt$l, header = T, sep = "\t", stringsAsFactors = F)
} else { stop(paste0("File '", opt$l, "' does not exist.")) }

# Reading the gene list file:
if(file.exists(opt$g)) {
    if(opt$v) print(paste0("Reading the gene list: ", opt$g))
    genes <- read.table(file = opt$g, header = F)
    colnames(genes) <- c("chr", "start", "stop", "name")
} else { stop(paste0("File '", opt$g, "' does not exist.")) }

# Checking if the top SNP is in the file:
if(opt$s %in% summ_stat$rsid) {
    rsid <- opt$s
    soi <- summ_stat[summ_stat$rsid == rsid,]
} else { stop(paste0("SNP: '", opt$s, "' is not in the summary file: '", opt$i ,"'"))}

# Checking the window size:
if(opt$w <= 0) {
    stop("The window size must be > 0.")
} else { win_size <- as.numeric(opt$w) }

min_pos <- min(soi$base_pair_location) - win_size/2
max_pos <- max(soi$base_pair_location) + win_size/2

# Merging LD info with summary statistics:
if(opt$v) print("Merging LD and summary statistics")
summ_stat_ld <- merge(x = summ_stat, 
                      y = LD_info[,c("SNP_B", "R2")],
                      by.x = "rsid", by.y = "SNP_B",
                      all.x = TRUE, all.y = FALSE)

summ_stat_sub <- summ_stat_ld[summ_stat_ld$chromosome == soi$chromosome & 
                              summ_stat_ld$base_pair_location > min_pos & 
                              summ_stat_ld$base_pair_location < max_pos,]

summ_stat_sub$R2[is.na(summ_stat_sub$R2)] <- 0

# Preparing gene data:
genes_in_region <- genes[genes$chr == soi$chromosome & 
                         genes$stop > min_pos & 
                         genes$start < max_pos,]
genes_in_region <- genes_in_region[order(genes_in_region$start),]
# Removing coordinates outside the range:
genes_in_region$start[genes_in_region$start < min_pos] <- min_pos
genes_in_region$stop[genes_in_region$stop > max_pos] <- max_pos 

# We reserve 20% of the y-axis to plot genes:
ymax <- max(-log10(summ_stat_sub$p_value))
ymin <- 0 - ymax/5
genes_y_pos <- rep(x = c(ymin, ymin/2, ymin/4), times = length(genes_in_region$name))[1:length(genes_in_region$name)]

# Plotting the manhattan plot around the top SNP:
gg_region <- ggplot(summ_stat_sub, aes(x = base_pair_location, y = -log10(p_value))) + 
               geom_point(data = soi, aes(x = base_pair_location, y = -log10(p_value)), colour = "purple", pch = 18, size = 5) +
               geom_text(data = soi, aes(x = base_pair_location, y = -log10(p_value), label = rsid), colour = "purple", hjust = -0.3, vjust = -0.1) +
               geom_segment(data = genes_in_region, aes(x = start, xend = stop, y = genes_y_pos, yend = genes_y_pos, colour = name), size = 4, alpha = 0.2) +
               geom_text(data = genes_in_region, aes(x = start, y = genes_y_pos, label = name, vjust = 1, colour = name)) +
               theme_bw() + scale_fill_gradient(low = "#414487", high = "#FDE725FF") + xlab(paste0("chr", soi$chromosome)) +
               ylim(c(ymin, ymax))
               
# Plot colored by R2:
gg_R2 <- gg_region + geom_point(aes(fill = R2), pch = 21, size = 3)

# Plot colored by Effect Size: can be "odds_ratio" or "beta" depending on the phenotype.
ES_name <- colnames(summ_stat_sub)[6]
if(ES_name == "odds_ratio") gg_ES <- gg_region + geom_point(aes(fill = odds_ratio), pch = 21, size = 3)
if(ES_name == "beta")       gg_ES <- gg_region + geom_point(aes(fill = beta), pch = 21, size = 3)
# (probably there is a more elegant way of dealing with a variable column name with ggplot2)

if(opt$v) print(paste0("Writing plot to: '", opt$o, "'.pdf"))
pdf(paste0(opt$o,".pdf"), width = 12)
    print(gg_R2)
    print(gg_ES)
dev.off()
