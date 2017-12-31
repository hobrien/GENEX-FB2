library(readr)
library(dplyr)
library(tools)
library(ggplot2)
library(optparse)

# Test for eGene enrichment (ratio of sig to non-sig GTEx eGenes that are sig in query sample vs. ratio in eGenes that are non-sig in query)
# I'm also determining the number of query topSNPs that are also sig in GTEx samples and the number that are sig for the same eGene
# It's not entirely clear (to me) what the appropriate background is to test for enrichment in these cases

option_list <- list(
  make_option(c("-o", "--outfile"), type="character", default="test.tsv", 
              help="output file"),
  make_option(c("-p", "--plot"), type="character", default="test.png", 
              help="plot of Discordance values")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser, positional_arguments=TRUE)
Discordance <- vector("list", length(opt$args))
for ( i in seq_along(opt$args) ) {
  GTcheckFile <- opt$args[[i]]
  sample<-file_path_sans_ext(basename(GTcheckFile))
  GTcheck <- read_tsv(GTcheckFile, 
                      col_names = c('CN', 'Discordance_total', 'Discordance_avg', 
                                    'Num_sites', 'SampleID', 'Sample_num'),
                      col_types = cols(
                        CN = col_character(),
                        Discordance_total = col_double(),
                        Discordance_avg = col_double(),
                        Num_sites = col_integer(),
                        SampleID = col_character(),
                        Sample_num = col_integer()
                      ),
                      comment = "#", trim_ws = TRUE)
  GTcheck <- GTcheck %>% mutate(refID=sample) %>% select(-CN, -Sample_num)
  Discordance[[i]] <- GTcheck
}

Discordance <- bind_rows(Discordance)
write_tsv(Discordance, opt$options$outfile)
ggplot(Discordance, aes(x=refID, y=Discordance_total)) + 
  geom_point() + 
  geom_point(data=filter(Discordance, SampleID==refID), colour='red')
ggsave(opt$options$plot)
