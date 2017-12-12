library(readr)
library(dplyr)
library(stringr)
library(optparse)

# Test for eGene enrichment (ratio of sig to non-sig GTEx eGenes that are sig in query sample vs. ratio in eGenes that are non-sig in query)
# I'm also determining the number of query topSNPs that are also sig in GTEx samples and the number that are sig for the same eGene
# It's not entirely clear (to me) what the appropriate background is to test for enrichment in these cases

option_list <- list(
  make_option(c("-q", "--query"), type="character",
              default="Shiny/GENEX-FB2/Data/results.bed", 
              help="query bed file"),
  make_option(c("-o", "--outfile"), type="character", default="test.txt", 
              help="output file")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser, positional_arguments=TRUE)
query <- read_tsv(opt$options$query, 
                  col_names=c("Chr", "start", "pos", "geneID", "cisVariants", "Beta1", 
                              "Beta2", "topSNP", "distance", "slope", "nominal_p", 
                              "nominal_p_threshold", "padj_direct", "padj_beta", "qvalue"))

overlaps <- data.frame()
for ( referenceFile in opt$args ) {
  sample <- str_replace(basename(referenceFile), '\\..*', '')
  referenceAll <- read_tsv(paste0('GTEx_Analysis_v7_eQTL/', sample, '.v7.egenes.txt.gz'))
  
  Combined <-  referenceAll %>% select(gene_id, gtex_qvalue=qval) %>%
    mutate(geneID=str_replace(gene_id, "\\.\\d+", "")) %>%
    inner_join(query)
  
  Overlap <- nrow(Combined)
  GTExSig <- filter(Combined, gtex_qvalue < .05) %>% nrow()
  QuerySig <- filter(Combined, qvalue<.05) %>% nrow()
  SigOverlap <- filter(Combined, qvalue<.05 & gtex_qvalue < .05) %>% nrow()
  test <- fisher.test(matrix(c(SigOverlap, 
                               GTExSig-SigOverlap,
                               QuerySig-SigOverlap,
                               Overlap-GTExSig-QuerySig+SigOverlap), nrow=2))
  
  referenceSig <- read_tsv(referenceFile, 
                           col_names=c("Chr", "start", "pos", "SNP", "geneID", "strand",
                                       "nominal_p", "slope", "slope_se", "qvalue"))
  
  SNPoverlap <- referenceSig %>% 
    select(Chr, pos, gtex_qvalue=qvalue) %>%
    inner_join(filter(query, qvalue<.05)) %>%
    group_by(Chr, pos) %>%
    slice(1) %>%
    nrow()
  
  SNP_eGene_overlap <-  referenceSig %>% select(Chr, pos, geneID, gtex_qvalue=qvalue) %>%
    mutate(geneID=str_replace(geneID, "\\.\\d+", "")) %>%
    inner_join(filter(query, qvalue<.05)) %>% nrow()
  
  overlaps <- overlaps %>% bind_rows(data.frame(sample=c(sample), 
                                                GTExSig=c(GTExSig), 
                                                QuerySig=c(QuerySig), 
                                                Overlap=c(Overlap), 
                                                SigOverlap=c(SigOverlap), 
                                                SNPoverlap=c(SNPoverlap),
                                                SNP_eGene_overlap= c(SNP_eGene_overlap),
                                                enrichment=c(test$estimate),
                                                p_enrich=c(test$p.value)
  )
  )
}

write_tsv(overlaps, opt$options$outfile)
