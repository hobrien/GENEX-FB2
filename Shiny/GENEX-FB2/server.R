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

# setwd("~/BTSync/FetalRNAseq/Github/GENEX-FB2/Shiny/GENEX-FB2/")

################################## Define functions ##################################

PlotEQTL<-function(row_num, counts, cis, target, snp_header) {
  qtl_stats <- cis[row_num,]
  geneID <- qtl_stats$geneID
  snp <- qtl_stats$topSNP
  snp_pos <- qtl_stats$pos
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
  #mean <- data %>% group_by(genotype) %>% summarise(mean = mean(value))
  title<-paste0(geneID, ' x ', snp)
  plot<-  ggplot(data, aes(y=value, x=DS)) + 
    #geom_errorbar(aes(ymin=mean, ymax=mean), colour='black', size=1, width=.5, data=mean) +
    geom_smooth(method='lm', colour='black') +
    geom_quasirandom(aes(colour=geno), width=0.05, alpha=.5) + 
    ylab("Normalised Counts") +
    xlab('Dosage') +
    #scale_x_discrete(breaks=c(0,1,2), labels=c('Ref', 'Het', 'Alt')) +
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

################################## Load Data ##################################

counts <- read_delim("./Data/counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  select(-Chr, -gene_type, -ChrType)

colnames(counts) <- str_replace_all(colnames(counts), 'norm.', '')

top_cis <- read_delim("./Data/cis_eqtl.txt", " ", escape_double = FALSE, trim_ws = TRUE, col_types = cols(pos='c')) %>%
  select(-one_of(c('cisVariants', 'Beta1', 'Beta2', 'nominal_p_threshold', 'padj_direct', 'padj_beta')))

target <- read_delim("./Data/SampleInfo.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(Sample=as.character(Sample))

snp_header <- system("bcftools view -h ./Data/combined_filtered.vcf.gz | tail -1", intern=TRUE) %>%
  str_split('\t')

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
  output$TopCisTable <- DT::renderDataTable({
    DT::datatable(top_cis[top_cis[, input$p_type] < input$pvalue, ], selection = 'single')
  })
  output$top_cis_download <- downloadHandler(
    filename = function() { 'cis_eqtl.csv' },
    content = function(file) {
      write_tsv(cis[cis$padj < input$pvalue, ], file)
    }
  )
  
  
})
