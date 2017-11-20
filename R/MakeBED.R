#setwd("~/BTSync/FetalRNAseq/Github/GENEX-FB2/")

library(tidyverse)
library(stringr)
library("optparse")

option_list <- list(
  make_option(c("--counts"), type="character", default="../Data/counts_vst.txt",
              help="Counts file from SARtools"),
  make_option(c("--genes"), type="character", default="../Data/geneloc.txt",
              help="Positions of genes (geneID, chr, start, end)"),
  make_option(c("--min"), type="numeric", default=6,
              help="Minimum count"),
  make_option(c("--num"), type="integer", default=10,
              help="Number of samples with minimum count"),
  make_option(c("--average"), action='store_true', type="logical", default=FALSE, 
              help="Use average expression rather than number of samples"),
  make_option(c("--out"),type="character", default="../Data/expression.bed", 
              help="Output file (BED)"),
  make_option(c("-e", "--exclude"), type="character", default='', 
              help="IDs of samples to exclude")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
exclude <- strsplit(opt$exclude, ',')[[1]]
# Import normalised counts
print(paste("reading", opt$counts))
counts <- read_delim(opt$counts, "\t", escape_double = FALSE, trim_ws = TRUE)
counts<-select(counts, -one_of(exclude))

if (opt$average){
  counts<-counts %>% filter(rowSums(.[,2:ncol(counts)]/(ncol(counts)-1))>opt$min)
} else{
  counts<-counts %>% filter(rowSums(.[,2:ncol(counts)] > opt$min)>opt$num)
}

print(paste("reading", opt$genes))
genepos <- read_delim(opt$genes, " ", col_names = c("id", "chr", "s1", "s2"), escape_double = FALSE, trim_ws = TRUE)
genepos <- genepos %>% mutate(id= str_replace(id, '(ENSG\\d+)\\.\\d+', '\\1')) %>% 
  inner_join(counts, by=c("id" = 'Id')) %>%
  filter(chr %in% paste0('chr', c(seq(22), 'X', 'Y'))) %>%
  select(`#Chr`=chr, start=s1, end=s2, ID=id, everything())

write_tsv(genepos, opt$out)