# Script to create the groupFile to input in SAIGE+ for the RVAT analysis:

# Takes the plink.annot.clean file created within the CreateGroupFile process and 
# writes a .groupFile for the RVAT analysis with SAIGE+

if(!require(optparse, quietly = T)) install.packages(optparse) 

option_list <- list(
    make_option(c("-a", "--annot"), 
                type = "character", 
                default = "", 
                help = "Full path to plink.annot.clean file.",
                metavar = "<PATH_TO_ANNOT>"),

    make_option(c("-v", "--verbose"),
                type = "logical",
                default = FALSE,
                help = "Print information during execution (FALSE by default)",
                metavar = "<TRUE,FALSE>")
);

opt <- parse_args(OptionParser(option_list = option_list))

# Reading the .fam file
if(file.exists(opt$a)) {
    if(opt$v) print(paste0("Reading the plink.annot.clean file", opt$a))
    annot <- read.table(file = opt$a, header = T, sep = ' ')
} else { stop(paste0("File '", opt$a, "' does not exist.")) }

#valid_annot <- annot[,2] != '.'
annot <- annot[annot$ANNOT != '.', ]
annot$anno <- 'no_annot'

colnames(annot) <- c("var", "ID", "anno") # We use the colnames to have the "correct" name for the 2nd column of the groupFile

annot_df <- as.data.frame(aggregate(ID ~ var + anno, annot, paste, collapse = ' '))

#groupFile <- as.data.frame(aggregate(ID ~ ANNOT, annot, paste, collapse = ' '))
# annot_df <- data.frame(group=c())
# for(gene in unique(annot$ANNOT)) {
#     annot_df <- rbind(annot_df, paste0(c(gene, 'var', annot$ID[annot$ANNOT == gene]), collapse = ' '))
#     annot_df <- rbind(annot_df, paste0(c(gene, 'anno', annot$IMPACT[annot$ANNOT == gene]), collapse = ' '))
# }

write.table(x = annot_df, file = 'saige.groupFile', quote = F, row.names = F, sep = ' ', col.names = F)
