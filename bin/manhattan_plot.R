# Script to create a manhattan plot from a SAIGE output
# Writes it to a pdf file

if(!require(optparse, quietly = T)) install.packages(optparse) 
if(!require(ggplot2, quietly = T))  install.packages(ggplot2) 
if(!require(qqman, quietly = T))    install.packages(qqman) 




qqman::manhattan(x = saige_sv, chr = "CHR", bp = "POS", snp = "MarkerID", p = "p.value", 
	 		     col = c("#6495ED", "#DE3163"), main = "Manhattan plot")
