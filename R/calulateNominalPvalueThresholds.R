suppressMessages(library(qvalue))

args <- commandArgs(trailingOnly = TRUE)

ifile = args[1]
fdr = as.numeric(args[2]);

cat("Processing fastQTL concatenated output [", ifile, "] controlling for FDR =", fdr * 100, "%\n");

#Read data
D = read.table(ifile, hea=FALSE, stringsAsFactors=FALSE)
colnames(D) <- c("geneID", "cisVariants", "Beta1", "Beta2", "Dummy", "topSNP", 
                 "distance", "nominal_p", "corrected_p_direct", "corrected_p_beta")
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
                 "distance", "nominal_p", "padj_direct", "padj_beta", "qval", "nominal_p_threshold")
#Write output
write.table(D[, c(1,2,3,4,6,7,8,12,9,10,11)], args[3], quote=FALSE, row.names=FALSE, col.names=TRUE)

cat("Done\n")
