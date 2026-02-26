#!/usr/bin/env Rscript

# Script to create the groupFile to input in SAIGE+ for the RVAT analysis:
# Takes the genes.annot.clean file created within the CreateGroupFile process and 
# writes a .groupFile for the RVAT analysis with SAIGE+

if(!require(optparse, quietly = T)) install.packages(optparse) 

option_list <- list(
    make_option(c("-g", "--genesets"), 
                type = "character", 
                default = "", 
                help = "Full path to 'genes.annot.clean' file.",
                metavar = "<PATH_TO_GENES>"),

    make_option(c("-a", "--annot"),
                type = "character",
                default = "",
                help = "Full path to the annotation file.",
                metavar = "<PATH_TO_ANNOT>"),

    make_option(c("-v", "--verbose"),
                type = "logical",
                default = FALSE,
                help = "Print information during execution (FALSE by default)",
                metavar = "<TRUE,FALSE>")
);

opt <- parse_args(OptionParser(option_list = option_list))

# Reading the .fam file
if(file.exists(opt$g)) {
    if(opt$v) print(paste0("Reading the gene file ", opt$g))
    genes <- read.table(file = opt$g, header = T, sep = ' ')
} else { stop(paste0("File '", opt$g, "' does not exist.")) }

genes <- genes[genes$ANNOT != '.', ]
genes$ann_default <- 'no_annot'
colnames(genes) <- c("var", "ID", "ann_default")

# If an annotation file was given in the config, we complete the annotations:
if(opt$a != "") {
    if(opt$v) print(paste0("Annotation file present, reading ", opt$a))

    annot <- read.table(file = opt$a, header = F, sep = " ", stringsAsFactor = FALSE)
    colnames(annot) <- c("var", "ann")

    genes_merged <- merge(x = genes, y = annot, by = "var", all.x = T, all.y = F)

    #Replacing NA annotations with "no_annot":
    genes_merged$ann[is.na(genes_merged$ann)] <- genes_merged$ann_default[is.na(genes_merged$ann)]
    genes <- genes_merged[,c("var", "ID", "ann")]
}

# Group by the gene name "ID":
genes_var <- aggregate(var ~ ID, data = genes, FUN = paste, collapse = " ")
genes_ann <- aggregate(ann ~ ID, data = genes, FUN = paste, collapse = " ")

genes_var$type <- "var"
genes_ann$type <- "anno"

colnames(genes_var) <- c("ID", "value", "type")
colnames(genes_ann) <- c("ID", "value", "type")

long_genes <- rbind(genes_var[,c(1,3,2)],
                    genes_ann[,c(1,3,2)])

# Reordering by gene names:
long_genes <- long_genes[order(long_genes$ID),]

write.table(x = long_genes, file = 'saige.groupFile', quote = F, row.names = F, sep = ' ', col.names = F)
