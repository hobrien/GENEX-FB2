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
cis <- read_delim("./Data/cis_eqtl.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
top_cis <- cis %>% group_by(SYMBOL) %>% arrange(pvalue) %>% slice(1) %>% ungroup()
trans <- read_delim("./Data/trans_eqtl.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
top_trans <- trans %>% group_by(SYMBOL, cisSYMBOL) %>% arrange(pvalue) %>% slice(1) %>% ungroup()

target <- read_delim("./Data/target.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(label=as.character(label))

PlotEQTL<-function(row_num, counts, cis, target, snps) {
  qtl_stats <- cis[row_num,]
  geneID <- qtl_stats$SYMBOL
  snp <- qtl_stats$snps
  data <- counts %>% filter(SYMBOL == geneID | Id == geneID) %>%  
    dplyr::select(-SYMBOL, -Id, -Chr) %>%
    gather() %>%
    separate(key, into=c('norm', 'label'), sep='[.]') %>%
    dplyr::select(label, value) %>%
    left_join(target)
  data<- snps %>% filter(id==snp) %>%
    dplyr::select(-id) %>%
    gather(label, genotype) %>%
    mutate(genotype = factor(genotype)) %>%
    inner_join(data)
  statistic <- qtl_stats$statistic[1]
  pval <- qtl_stats$pvalue[1]
  qval <- qtl_stats$padj[1]
  mean <- data %>% group_by(genotype) %>% summarise(mean = mean(value))
  title<-paste0(geneID, ' x ', snp)
  plot<-  ggplot(data, aes(x=genotype, colour=genotype)) + 
    geom_errorbar(aes(ymin=mean, ymax=mean), colour='black', size=1, width=.5, data=mean) +
    geom_jitter(aes(y=value), height = 0, width=.1, alpha=.75) + 
    ylab("normalised counts") +
    xlab('') +
    scale_x_discrete(breaks=c(0,1,2), labels=c('Ref', 'Het', 'Alt')) +
    side_theme() +
    scale_colour_brewer(type = "qual", palette = 6) +
  ggtitle(title) 
  plot
}


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  output$SciNotation <- reactive({
    paste0("p-value: ", 10^input$pvalue)
  })
  output$eQTLplotTop <- renderPlot({
    PlotEQTL(input$TopCisTable_rows_selected, counts, top_cis[top_cis[, input$p_type] < 10^input$pvalue, ], target, snps)
  })
  output$eQTLplotAll <- renderPlot({
    PlotEQTL(input$AllCisTable_rows_selected, counts, cis[cis[, input$p_type] < 10^input$pvalue, ], target, snps)
  })
  output$eQTLplotTransTop <- renderPlot({
    PlotEQTL(input$TopTransTable_rows_selected, counts, top_trans[top_trans[, input$p_type] < 10^input$pvalue, ], target, snps)
  })
  output$eQTLplotTransAll <- renderPlot({
    PlotEQTL(input$AllTransTable_rows_selected, counts, trans[trans[, input$p_type] < 10^input$pvalue, ], target, snps)
  })
  output$TopCisTable <- DT::renderDataTable({
    DT::datatable(top_cis[top_cis[, input$p_type] < 10^input$pvalue, ], selection = 'single')
  })
  output$AllCisTable <- DT::renderDataTable({
    DT::datatable(cis[cis[, input$p_type] < 10^input$pvalue, ], selection = 'single')
  })
  output$TopTransTable <- DT::renderDataTable({
    DT::datatable(top_trans[top_trans[, input$p_type] < 10^input$pvalue, ], selection = 'single')
  })
  output$AllTransTable <- DT::renderDataTable({
    DT::datatable(trans[trans[, input$p_type] < 10^input$pvalue, ], selection = 'single')
  })
  output$download12_19 <- downloadHandler(
    filename = function() { 'cis_eqtl.csv' },
    content = function(file) {
      write_tsv(cis[cis$padj < input$pvalue, ], file)
    }
  )
  
  
})
