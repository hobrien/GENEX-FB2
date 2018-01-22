library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(argparser)

p <- arg_parser("Combines SNP info, gene info and eQTL info for SMR")
p <- add_argument(p, "gwas_summary", help="")
p <- add_argument(p, "outfile", help="")

args <- parse_args(p)

gwas_summary_file <- args$gwas_summary
out_file <- args$outfile

CLOZUK <- read_tsv(gwas_summary_file,
                   col_names=TRUE, cols(
                     SNP = col_character(),
                     Freq.A1 = col_double(),
                     CHR = col_integer(),
                     BP = col_integer(),
                     A1 = col_character(),
                     A2 = col_character(),
                     OR = col_double(),
                     SE = col_double(),
                     P = col_double(),
                     Ntot = col_integer()
                   ), trim_ws = TRUE)



CLOZUK <- separate(CLOZUK, SNP, c('chr', 'pos', 'A1.1', 'A2.2'), extra='merge') %>%
  dplyr::filter(CLOZUK, str_detect(chr, '^rs')) %>%
  filter(!chr %in% CLOZUK_rsID$chr[duplicated(CLOZUK_rsID$chr)]) %>%
  select(SNP=chr, A1, A2, freq=Freq.A1, b=OR, se=SE, p=P, n=Ntot) %>%
  mutate(A1=toupper(A1), A2=toupper(A2))

write_tsv(CLOZUK, out_file, col_names=TRUE)
