#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(purrr)
source("FormatGGplot.R")
library(DT)
library(ggbeeswarm)

################################## Define functions ##################################
PlotEQTL<-function(row_num, counts, cis, target, snp_header) {
  qtl_stats <- cis[row_num,]
  geneID <- qtl_stats$Id
  snp <- qtl_stats$topSNP
  snp_pos <- paste(qtl_stats$Chr, qtl_stats$pos, sep=":")
  expression <- counts %>% filter(Id == geneID) %>%  
    dplyr::select(-Id, -SYMBOL) %>%
    gather("Sample", "value") %>%
    dplyr::select(Sample, value) %>%
    left_join(target)
  genotypes <- system(paste("bcftools view -H -r", snp_pos, "./Data/combined_filtered.vcf.gz"), intern=TRUE) %>%
    str_split('\t')
  ref <- genotypes[[1]][4]
  alt <- genotypes[[1]][5]
  data <- tibble(Sample=snp_header[[1]][10:length(genotypes[[1]])], geno=genotypes[[1]][10:length(genotypes[[1]])]) %>%
    separate(geno, c("GT", "DS", "GP" ), sep=":") %>%
    full_join(expression) %>%
    filter(!is.na(GT)) %>%
    mutate(DS=as.numeric(DS), 
           geno=factor(ifelse(GT=='0|0', paste0(UQ(ref), UQ(ref)),
                              ifelse(GT=='1|1', paste0(UQ(alt), UQ(alt)),
                                     paste0(UQ(ref), UQ(alt)))),
                       levels=c(paste0(UQ(ref), UQ(ref)), paste0(UQ(ref), UQ(alt)), paste0(UQ(alt), UQ(alt)) )))
  pval <- qtl_stats$nominal_p[1]
  qval <- qtl_stats$qvalue[1]
  title<-paste0(geneID, ' x ', snp)
  plot<-  ggplot(data, aes(y=value, x=DS)) + 
    geom_smooth(method='lm', colour='black') +
    geom_quasirandom(aes(colour=geno), width=0.05, alpha=.5) + 
    ylab("Normalised Counts") +
    xlab('Dosage') +
    side_theme() +
    scale_colour_brewer(type = "qual", palette = 6) +
    ggtitle(title) +
    theme(legend.position="right") +
    labs(color='Best Guess\nGenotype')
  plot
}

add_links <-function(fitted) {
  mutate(fitted, SYMBOL=paste0("<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=", SYMBOL, " target='_blank'>", SYMBOL, "</a>"),
         GTEx=paste0("<a href=https://gtexportal.org/home/gene/", Id, " target='_blank'>GTEx</a>"),
         Id=paste0("<a href=http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=", Id, " target='_blank'>", Id, "</a>")
  )         
}

filter_table <- function(fitted, p_type, p_val) {
  fitted <- filter(fitted, !is.na(UQ(as.name(p_type))) & UQ(as.name(p_type)) <= p_val ) 
  fitted <- mutate(fitted, nominal_p = as.numeric(format(nominal_p, digits=3)), qvalue = as.numeric(format(qvalue, digits=3))) %>%
    dplyr::select(Id, SYMBOL, `Top SNP`=topSNP, Chr, `SNP Pos`=pos, slope, `pvalue`=nominal_p, padj=qvalue)
  fitted
}

################################## Load Data ##################################
target <- read_delim("./Data/SampleInfo.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(Sample=as.character(Sample))

snp_header <- system("bcftools view -h ./Data/combined_filtered.vcf.gz | tail -1", intern=TRUE) %>%
  str_split('\t')

# Gene-level analysis
counts <- read_delim("./Data/counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(-Chr, -gene_type, -ChrType)

colnames(counts) <- str_replace_all(colnames(counts), 'norm.', '')

top_cis <- read_delim("./Data/egenes_gene_q05.bed.gz", "\t", escape_double = FALSE, trim_ws = TRUE, 
                      col_names=c("Chr", "start", "pos", "geneID", "cisVariants", "Beta1", 
                                  "Beta2", "true_df", "pval_true_df", "topSNP", "distance", 
                                  " minor_allele_samples", "minor_allele_count", "maf", 
                                  "ref_factor", "nominal_p", "slope", "slope_se", "padj_direct", 
                                  "padj_beta", "qvalue", "nominal_p_threshold")) %>%
  left_join(dplyr::select(counts, Id, SYMBOL), by=c('geneID' = 'Id')) %>%
  dplyr::select(-one_of(c('start', 'cisVariants', 'Beta1', 'Beta2', 'nominal_p_threshold', 'nominal_p_threshold.1','padj_direct', 'padj_beta'))) %>%
  mutate(Chr=str_replace(Chr, 'chr', '')) %>%
  dplyr::rename(Id=geneID)

# Transcript-level analysis
counts_tr <- read_delim("./Data/counts_tr.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(-one_of(c("Chr", "gene_type", "ChrType", "gene_id")))

colnames(counts_tr) <- str_replace_all(colnames(counts_tr), 'norm.', '')

top_cis_tr <- read_delim("./Data/egenes_transcript_q05.bed.gz", "\t", escape_double = FALSE, trim_ws = TRUE, 
                      col_names=c("Chr", "start", "pos", "geneID", "cisVariants", "Beta1", 
                                  "Beta2", "true_df", "pval_true_df", "topSNP", "distance", 
                                  " minor_allele_samples", "minor_allele_count", "maf", 
                                  "ref_factor", "nominal_p", "slope", "slope_se", "padj_direct", 
                                  "padj_beta", "qvalue", "nominal_p_threshold")) %>%
  dplyr::rename(Id=geneID) %>%
  left_join(dplyr::select(counts_tr, Id, SYMBOL), by=c('Id')) %>%
  dplyr::select(-one_of(c('geneID', 'start', 'cisVariants', 'Beta1', 'Beta2', 'nominal_p_threshold', 'nominal_p_threshold.1','padj_direct', 'padj_beta'))) %>%
  mutate(Chr=str_replace(Chr, 'chr', '')) 
  

################################## Run server ##################################

shinyServer(function(session, input, output) {
  observe({
    updateSliderInput(session, "pvalue", value = input$typedPval)
  })
  output$eQTLplotTop <- renderPlot({
    validate(
      need(input$TopCisTable_rows_selected != "", "Please select a row from the table")
    )
    PlotEQTL(input$TopCisTable_rows_selected, counts, top_cis[top_cis[, input$p_type] < input$pvalue, ], target, snp_header)
  })
  output$eQTLplotTopTr <- renderPlot({
    validate(
      need(input$TopCisTableTr_rows_selected != "", "Please select a row from the table")
    )
    PlotEQTL(input$TopCisTableTr_rows_selected, counts_tr, top_cis_tr[top_cis_tr[, input$p_type] < input$pvalue, ], target, snp_header)
  })
  output$TopCisTable <- DT::renderDataTable({
    DT::datatable(add_links(filter_table(top_cis, input$p_type, input$pvalue)), escape = FALSE, selection="single", caption = 'Top eQTL for each eGENE')
  })
  output$TopCisTableTr <- DT::renderDataTable({
    DT::datatable(add_links(filter_table(top_cis_tr, input$p_type, input$pvalue)), escape = FALSE, selection="single", caption = 'Top eQTL for each eTRANSCRIPT')
  })
 })

