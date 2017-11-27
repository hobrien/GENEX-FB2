suppressMessages(library(qvalue))
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

ifile = args[1]
snp_file = args[2]
fdr = as.numeric(args[3]);

cat("Processing fastQTL concatenated output [", ifile, "] controlling for FDR =", fdr * 100, "%\n");

#Read data
D = read.table(ifile, hea=FALSE, stringsAsFactors=FALSE)

D = D[which(!is.na(D[, 10])),]
cat("  * Number of molecular phenotypes =" , nrow(D), "\n")
cat("  * Correlation between Beta approx. and Empirical p-values =", round(cor(D[, 9], D[, 10]), 4), "\n")

#Run qvalue on pvalues for best signals
Q = qvalue(D[, 10])
D$qval = Q$qvalue
cat("  * Proportion of significant phenotypes =" , round((1 - Q$pi0) * 100, 2), "%\n")

#Determine significance threshold
set0 = D[which(D$qval <= fdr),] 
set1 = D[which(D$qval > fdr),]
pthreshold = (sort(set1$V10)[1] - sort(-1.0 * set0$V10)[1]) / 2
cat("  * Corrected p-value threshold = ", pthreshold, "\n")

#Calculate nominal pvalue thresholds
D$nthresholds = qbeta(pthreshold, D$V3, D$V4, ncp = 0, lower.tail = TRUE, log.p = FALSE)

colnames(D) <- c("geneID", "cisVariants", "Beta1", "Beta2", "Dummy", "topSNP", 
                 "distance", "nominal_p", "slope", "padj_direct", "padj_beta", "qvalue", "nominal_p_threshold")
# Add SNP positions
snp_pos <- read_tsv(args[2], col_names=FALSE)

snp_pos <- mutate(snp_pos, pos = paste(X1, X2, sep=':')) %>%
  select(topSNP=X3, pos)
D <- left_join(D, snp_pos)

#Write output
write.table(D[, c(1,2,3,4,6,7,9,8,13,10,11,12,14)], args[4], quote=FALSE, row.names=FALSE, col.names=TRUE)

cat("Done\n")
