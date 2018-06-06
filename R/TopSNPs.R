suppressMessages(library(dplyr))
suppressMessages(library(readr))

suppressMessages(library(argparser))

# parse inputs
p <- arg_parser("Selects the lowest p-value for each SNP")
p <- add_argument(p, "input", help="")
p <- add_argument(p, "snpfile", help="")
p <- add_argument(p, "output", help="")
args <- parse_args(p)

all_snps <- read_tsv(args$input, col_names=c("gene_id", "variant_id", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se"))
sig_snps <- all_snps %>% arrange(pval_nominal) %>% 
  group_by(variant_id) %>% 
  slice(1) %>% 
  ungroup()

snp_pos <- read_tsv(args$snpfile, col_names=FALSE)

snp_pos <- snp_pos %>% select(variant_id=X3, A1=X4, A2=X5)
sig_snps <- left_join(sig_snps, snp_pos) %>%
  select(SNP=variant_id, A1, A2, freq=maf, b=slope, se=slope_se, p=pval_nominal)
write_tsv(sig_snps, args$output)
