library(readr)
library(dplyr)
library(stringr)
library(argparser)

p <- arg_parser("Combines SNP info, gene info and eQTL info for SMR")
p <- add_argument(p, "fastqtlOutput", help="")
p <- add_argument(p, "snpfile", help="")
p <- add_argument(p, "genesfile", help="")
p <- add_argument(p, "outfile", help="")

args <- parse_args(p)

eqtl_file <- args$fastqtlOutput
snp_file <- args$snpfile
gene_file <- args$genesfile
out_file <- args$outfile

snp_pos <- read_tsv(snp_file, col_names = c('Chr', 'BP', 'SNP', 'A1', 'A2'), # need to add A1, A2
                       trim_ws = TRUE)

snp_pos <- mutate(snp_pos, chr=paste0('chr', Chr))

gene_pos <- read_tsv(gene_file, col_names = c('Gene', 'Probe_Chr', 'start', 'end', 'Orientation'),
                    trim_ws = TRUE)
gene_pos <- gene_pos %>% mutate(Probe_bp = ifelse(Orientation == '+', start, end), Gene = str_replace(Gene, '\\.\\d+', ''))

eqtls <- read_tsv(eqtl_file, col_names = c('Gene', 'SNP', 'tss_distance', 'ma_samples', 
                                           'ma_count', 'Freq', 'p', 'b', 'se'), trim_ws=TRUE)

combined <- eqtls %>% left_join(snp_pos) %>%
  left_join(gene_pos) %>%
  select(SNP, Chr, BP, A1, A2, Freq, Probe=Gene, Probe_Chr, Probe_bp, Gene, Orientation, b, se, p)

write_tsv(combined, out_file)
