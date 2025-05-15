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

# Group by the gene name "ID":
annot <- aggregate(x = annot, by = list(annot$ID, annot$anno), FUN = paste, simplify = F)
# Then group the dataframe by genes to have the variants and annotations and two different lines:
annot <- annot[,c("Group.1", "var", "anno")]

# Transforming lists into characters:
annot$var <- as.character(annot$var)
annot$anno <- as.character(annot$anno)

# Removing unwanted characters (quotes and "c()"):
annot$var <- paste(gsub(pattern = '"',    replacement = '', x = annot$var))
annot$var <- paste(gsub(pattern = ',',    replacement = '', x = annot$var))
annot$var <- paste(gsub(pattern = '\n',   replacement = '', x = annot$var))
annot$var <- paste(gsub(pattern = 'c\\(', replacement = '', x = annot$var))
annot$var <- paste(gsub(pattern = '\\)$', replacement = '', x = annot$var))

annot$anno <- paste(gsub(pattern = '"',    replacement = '', x = annot$anno))
annot$anno <- paste(gsub(pattern = ',',    replacement = '', x = annot$anno))
annot$anno <- paste(gsub(pattern = '\n',   replacement = '', x = annot$anno))
annot$anno <- paste(gsub(pattern = 'c\\(', replacement = '', x = annot$anno))
annot$anno <- paste(gsub(pattern = '\\)$', replacement = '', x = annot$anno))

# Transforming the annotations from wide to long format:
long_annot <- reshape(data = annot, 
                      direction = "long", 
                      v.names = "value", 
                      varying = c("var", "anno"), 
                      timevar = "type", 
                      times = names(annot)[2:3])
long_annot <- long_annot[,c(1,2,3)]

# Reordering by gene names:
long_annot <- long_annot[order(long_annot$Group.1),]

write.table(x = long_annot, file = 'saige.groupFile', quote = F, row.names = F, sep = ' ', col.names = F)
