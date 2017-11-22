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
source("FormatGGplot.R")
library(DT)



# setwd("~/BTSync/FetalRNAseq/Github/GENEX-FB2/Shiny/GENEX-FB2/")
counts <- read_delim("./Data/counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
top_cis <- read_delim("./Data/cis_eqtl.txt", " ", escape_double = FALSE, trim_ws = TRUE)

target <- read_delim("./Data/SampleInfo.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(Sample=as.character(Sample))

snp_header <- system("bcftools view -h ./Data/genotypes.vcf.gz | tail -1", intern=TRUE) %>%
  str_split('\t')

PlotEQTL<-function(row_num, counts, cis, target, snp_header) {
  qtl_stats <- cis[row_num,]
  geneID <- qtl_stats$geneID
  snp <- qtl_stats$topSNP
  snp_pos <- qtl_stats$snp_pos
  data <- counts %>% filter(Id == geneID) %>%  
    dplyr::select(-Id) %>%
    gather("Sample", "value") %>%
    dplyr::select(Sample, value) %>%
    left_join(target)
  genotypes <- system(paste("bcftools view -H -r", snp_pos, "./Data/genotypes.vcf.gz"), intern=TRUE) %>%
    str_split('\t') %>%
    map(str_extract())
  data <- as.tibble(matrix(c(snp_header, genotypes), nrow=2)) %>%
    select(-pos, -snp) %>%
    gather(SampleID, dose) %>%
  statistic <- qtl_stats$statistic[1]
  pval <- qtl_stats$nominal_p[1]
  qval <- qtl_stats$qvalue[1]
  #mean <- data %>% group_by(genotype) %>% summarise(mean = mean(value))
  title<-paste0(geneID, ' x ', snp)
  plot<-  ggplot(data, aes(x=genotype, colour=genotype)) + 
    #geom_errorbar(aes(ymin=mean, ymax=mean), colour='black', size=1, width=.5, data=mean) +
    geom_jitter(aes(y=value), height = 0, width=.1, alpha=.75) + 
    ylab("normalised counts") +
    xlab('') +
    scale_x_discrete(breaks=c(0,1,2), labels=c('Ref', 'Het', 'Alt')) +
    side_theme() +
    scale_colour_brewer(type = "qual", palette = 6) +
  ggtitle(title) 
  plot
}

shinyServer(function(input, output) {
  output$SciNotation <- reactive({
    paste0("p-value: ", 10^input$pvalue)
  })
  output$eQTLplotTop <- renderPlot({
    PlotEQTL(input$TopCisTable_rows_selected, counts, top_cis[top_cis[, input$p_type] < 10^input$pvalue, ], target, snps)
  })
  output$TopCisTable <- DT::renderDataTable({
    DT::datatable(top_cis[top_cis[, input$p_type] < 10^input$pvalue, ], selection = 'single')
  })
  output$top_cis_download <- downloadHandler(
    filename = function() { 'cis_eqtl.csv' },
    content = function(file) {
      write_tsv(cis[cis$padj < input$pvalue, ], file)
    }
  )
  
  
})
