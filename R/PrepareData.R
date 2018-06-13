library(readr)
library(dplyr)
library("optparse")

option_list <- list (
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="Name of output file"),
  make_option(c("-l", "--lambda"), type="numeric", default=NULL, 
              help=""),
  make_option(c("-s", "--sample_info"), type="character", default=NULL, 
              help=""),
  make_option(c("-v", "--vcf"), type="character", default=NULL, 
              help=""),
  make_option(c("-l", "--lambda"), type="numeric", default=NULL, 
              help="")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser, positional_arguments=FALSE)

gene_info <- read_tsv("../Data/genes.txt") %>%
  mutate(gene_id = sub("\\.[0-9]+", "", gene_id)) %>%
  dplyr::select(Id = gene_id, SYMBOL=gene_name, Chr=seqid)

cis <- read_delim("../MatrixEQTL/cis_eqt.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
cis <- mutate(cis, gene = sub("(ENSG[0-9]+)\\.[0-9]+", '\\1', gene),
         statistic = as.numeric(format(statistic, digits=2)), 
         pvalue = as.numeric(format(pvalue, digits=2)), 
         FDR = as.numeric(format(FDR, digits=2)),
         beta = as.numeric(format(beta, digits=2))) %>%
  rename(Id=gene, padj=FDR)

cis <- right_join(gene_info, cis)
write_tsv(cis, "../Shiny/GENEX-FB2/Data/cis_eqtl.txt")

trans <- read_delim("../MatrixEQTL/trans_eqt.txt", "\t", escape_double = FALSE, trim_ws = TRUE) 
trans <- mutate(trans, gene = sub("(ENSG[0-9]+)\\.[0-9]+", '\\1', gene),
         statistic = as.numeric(format(statistic, digits=2)), 
         pvalue = as.numeric(format(pvalue, digits=2)), 
         FDR = as.numeric(format(FDR, digits=2)),
         beta = as.numeric(format(beta, digits=2))) %>%
  rename(Id=gene, padj=FDR)

trans <- dplyr::select(cis, snps, cisId=Id, cisSYMBOL=SYMBOL) %>% left_join(trans)
trans <- right_join(gene_info, trans)
write_tsv(trans, "../Shiny/GENEX-FB2/Data/trans_eqtl.txt")

genotypes <- read_delim("../MatrixEQTL/genotypes_formatted.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
genotypes_filtered <- semi_join(genotypes, cis, by=c("id" = 'snps'))
write_tsv(genotypes_filtered, "../Shiny/GENEX-FB2/Data/genotypes.txt")

file.copy(opt$options$sample_info, "../Shiny/GENEX-FB2/Data/SampleInfo.txt", overwrite=TRUE)
file.copy(opt$options$vcf, "../Shiny/GENEX-FB2/Data/combined_filtered.vcf.gz", overwrite=TRUE)

file.copy("../Data/MalevsFemale.complete.txt", "../Shiny/GENEX-FB2/Data/counts.txt", overwrite=TRUE)
file.copy("../Data/fitted.txt", "../Shiny/GENEX-FB2/Data/fitted.txt", overwrite=TRUE)
