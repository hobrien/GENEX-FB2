suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(qvalue))
suppressMessages(library(tools))
library("optparse")

option_list <- list (
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="Name of output file"),
  make_option(c("-f", "--fdr"), type="numeric", default=0.05, 
              help=""),
  make_option(c("-s", "--snps"), type="character", default=NULL, 
              help=""),
  make_option(c("-l", "--lambda"), type="numeric", default=NULL, 
              help="")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser, positional_arguments=TRUE)


fastqtl_results_list=vector("list", length(opt$args))
for (i in seq_along(fastqtl_results_list)) {
  cat("Processing FastQTL output (", opt$args[i], "), with FDR=", opt$options$fdr, "\n", sep="")
  
  # input files have no headers
  D <- read.table(opt$args[i], header=FALSE, stringsAsFactors=FALSE)
  if (dim(D)[2]==17) {
    colnames(D) <- c('gene_id', 'num_var', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df', 'variant_id', 'tss_distance',
                     'minor_allele_samples', 'minor_allele_count', 'maf', 'ref_factor',
                     'pval_nominal', 'slope', 'slope_se', 'pval_perm', 'pval_beta')
  } else {
    stop("FastQTL output in unrecognized format (mismatched number of columns).")
  }
  fastqtl_results_list[[i]] <- D
}

D <- bind_rows(fastqtl_results_list)

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
if (is.null(opt$options$lambda) || is.na(opt$options$lambda)) {
  Q <- qvalue(D[, 'pval_beta'])
} else {
  cat("  * Calculating q-values with lambda = ", opt$options$lambda, "\n", sep="")
  Q <- qvalue(D[, 'pval_beta'], lambda=opt$options$lambda)
}

D$qval <- signif(Q$qvalues, 6)
cat("  * Proportion of significant phenotypes (1-pi0): " , round((1 - Q$pi0), 2), "\n", sep="")
cat("  * eGenes @ FDR ", opt$options$fdr, ":   ", sum(D[, 'qval']<opt$options$fdr), "\n", sep="")

# determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
ub <- sort(D[D$qval > opt$options$fdr, 'pval_beta'])[1]  # smallest p-value above FDR
lb <- -sort(-D[D$qval <= opt$options$fdr, 'pval_beta'])[1]  # largest p-value below FDR
pthreshold <- (lb+ub)/2
cat("smallest p-value above FDR: ", ub, "\n")
cat("largest p-value below FDR: ", lb, "\n")
cat("  * min p-value threshold @ FDR ", opt$options$fdr, ": ", pthreshold, "\n", sep="")
D[, 'pval_nominal_threshold'] <- signif(qbeta(pthreshold, D[, 'beta_shape1'], D[, 'beta_shape2'], ncp=0, lower.tail=TRUE, log.p=FALSE), 6)

snp_pos <- read_tsv(opt$options$snps, col_names=c('chr', 'start', 'end', 'variant_id', 'score', 'strand', 'A1', 'A2'))

snp_pos <- snp_pos %>% select(variant_id, chr, start, end)
D <- left_join(D, snp_pos) %>%
  select(chr, start, end, everything())

sig_egenes <- filter(D, qval<=opt$options$fdr)

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

write.table(sig_egenes, gzfile(opt$options$out), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
