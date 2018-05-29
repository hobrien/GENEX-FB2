#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyBS)

# Application title
#titlePanel("Gene Expression in the Fetal Brain: Sex Biases"),

navbarPage("Fetal Brain Sequencing (FBSeq) 1: eQTLs",
            tabPanel("Cis eQTLs",
                     sidebarLayout(
                       sidebarPanel(
                         radioButtons("p_type", "Maximum p-value", c('Uncorrected p-values' = 'nominal_p', 'FDR corrected p-values (q-values)'= 'qvalue'), selected = 'qvalue', inline = FALSE,
                                      width = NULL),
                         sliderInput("pvalue", "p-value:", 
                                     min = 0, max = 1, value = 0.1, step= 0.01),
                         textInput("typedPval", "Type p-value", value=.1),
                         conditionalPanel(
                           'input.dataset === "Gene-level analysis"',
                           HTML("<strong>Select row to plot data</strong><br>"),
                           actionButton("PlotTopCis", "Plot eQTL")
                         ),
                         conditionalPanel(
                           'input.dataset === "Transcript-level analysis"',
                           HTML("<strong>Select row to plot data</strong><br>"),
                           actionButton("PlotTopCisTr", "Plot eQTL")
                         )
                        ),
                       mainPanel(
                         tabsetPanel(
                           id = 'dataset',
                           tabPanel('Gene-level analysis', 
                                    DT::dataTableOutput('TopCisTable')),
                           bsModal("TopCisPlot", "Top Cis eQTL", "PlotTopCis", size = "large",plotOutput("eQTLplotTop")),
                           tabPanel('Transcript-level analysis', 
                                    DT::dataTableOutput('TopCisTableTr')),
                           bsModal("TopCisPlotTr", "Top Cis eQTL", "PlotTopCisTr", size = "large",plotOutput("eQTLplotTopTr"))
                         )   
                       )
                       
                     )
            )
)
