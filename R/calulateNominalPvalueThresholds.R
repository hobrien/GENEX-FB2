suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(qvalue))
suppressMessages(library(tools))
suppressMessages(library(argparser))

# parse inputs
p <- arg_parser("Annotates FastQTL permutation output and runs qvalue")
p <- add_argument(p, "fastqtlOutput", help="")
p <- add_argument(p, "snpfile", help="")
p <- add_argument(p, "fdr", type="numeric", help="")
p <- add_argument(p, "outfile", help="")
p <- add_argument(p, "filtered", help="")
#p <- add_argument(p, "lambda", type="numeric", help="", default=NULL)
args <- parse_args(p)

cat("Processing FastQTL output (", args$fastqtlOutput, "), with FDR=", args$fdr, "\n", sep="")

# input files have no headers
D <- read.table(args$fastqtlOutput, header=FALSE, stringsAsFactors=FALSE)
if (dim(D)[2]==17) {
  colnames(D) <- c('gene_id', 'num_var', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df', 'variant_id', 'tss_distance',
                   'minor_allele_samples', 'minor_allele_count', 'maf', 'ref_factor',
                   'pval_nominal', 'slope', 'slope_se', 'pval_perm', 'pval_beta')
} else {
  stop("FastQTL output in unrecognized format (mismatched number of columns).")
}

# remove duplicates (
# keep eGene with lowest nominal p-value (only applies when top_snp not tested in one dup)
# if nominal p identical (because top SNP tested in both), keep lowest corrected p-value (because window includes more SNPs)
D<-arrange(D, pval_nominal, desc(pval_beta)) %>% group_by(gene_id) %>% dplyr::slice(1) %>% ungroup() %>% as.data.frame()

# remove genes w/o variants
nanrows <- is.na(D[, 'pval_beta'])
D <- D[!nanrows, ]
cat("  * Number of genes tested: ", nrow(D), " (excluding ", sum(nanrows), " genes w/o variants)\n", sep="")
cat("  * Correlation between Beta-approximated and empirical p-values: ", round(cor(D[, 'pval_perm'], D[, 'pval_beta']), 4), "\n", sep="")

# calculate q-values
#if (is.null(args$lambda) || is.na(args$lambda)) {
  Q <- qvalue(D[, 'pval_beta'])
#} else {
#  cat("  * Calculating q-values with lambda = ", args$lambda, "\n", sep="")
#  Q <- qvalue(D[, 'pval_beta'], lambda=args$lambda)
#}

D$qval <- signif(Q$qvalues, 6)
cat("  * Proportion of significant phenotypes (1-pi0): " , round((1 - Q$pi0), 2), "\n", sep="")
cat("  * eGenes @ FDR ", args$fdr, ":   ", sum(D[, 'qval']<args$fdr), "\n", sep="")

# determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
ub <- sort(D[D$qval > args$fdr, 'pval_beta'])[1]  # smallest p-value above FDR
lb <- -sort(-D[D$qval <= args$fdr, 'pval_beta'])[1]  # largest p-value below FDR
pthreshold <- (lb+ub)/2
cat("smallest p-value above FDR: ", ub, "\n")
cat("largest p-value below FDR: ", lb, "\n")
cat("  * min p-value threshold @ FDR ", args$fdr, ": ", pthreshold, "\n", sep="")
D[, 'pval_nominal_threshold'] <- signif(qbeta(pthreshold, D[, 'beta_shape1'], D[, 'beta_shape2'], ncp=0, lower.tail=TRUE, log.p=FALSE), 6)

snp_pos <- read_tsv(args$snpfile, col_names=FALSE)

snp_pos <- mutate(snp_pos, chr=paste0('chr', X1), start = X2-1) %>%
  select(variant_id=X3, chr, start, end=X2)
D <- left_join(D, snp_pos) %>%
  select(chr, start, end, everything())

write.table(D, gzfile(args$outfile), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

sig_egenes <- filter(D, qval<=args$fdr)

# pval_nominal_threshold refers to pval_true_df, but the nominal pass of fastQTL only reports pval_nominal
# in most cases, pval_true_df is larger than pval_nominal, so filtering using 
# pval_nominal <= pval_nominal_threshold will result in non-sig SNPs being included (this
# is mitigated by only including eQTLs from sig eGenes)
# HOWEVER, occasionally pval_nominal is larger than pval_true_df, meaning that sig eQTLs
# are filtered out. There's not a lot I can do about this short calculating pval_true_df
# for everything, but I can modify pval_nominal_threshold for these genes so at least the
# top eQTL is not filtered. Otherwise, the number of eGenes in sig_eqtls will be less than
# the number in the eGenes file.
# Strangely, there doesn't appear to be any relationship between true_df and pval_true_df,
# so I have no idea how this is calculated

if(nrow(filter(sig_egenes, pval_nominal > pval_true_df)) >0) {
    cat(nrow(filter(sig_egenes, pval_nominal > pval_true_df)), " eGenes with df_corrected p-value < nominal p-value\n", sep="")
    cat(nrow(filter(sig_egenes, pval_nominal > pval_nominal_threshold)), " sig eGenes have nominal p > threshold. Modifying threshold so top eQTL will be sig\n", sep="") 
    sig_egenes <- sig_egenes %>% mutate(pval_nominal_threshold=ifelse(pval_nominal_threshold >= pval_nominal, pval_nominal_threshold, pval_nominal))
}

write.table(sig_egenes, gzfile(args$filtered), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
