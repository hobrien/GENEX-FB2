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
                           actionButton("PlotToCis", "Plot")
                         ),
                         conditionalPanel(
                           'input.dataset === "Transcript-level analysis"',
                           HTML("<strong>Select row to plot data</strong><br>"),
                           actionButton("PlotToCisTr", "Plot")
                         )
                        ),
                       mainPanel(
                         tabsetPanel(
                           id = 'dataset',
                           tabPanel('Gene-level analysis', 
                                    DT::dataTableOutput('TopCisTable'),
                                    DT::dataTableOutput('TransTable2')),
                           bsModal("ToCisPlot", "Top Cis eQTL", "PlotToCis", size = "large",plotOutput("eQTLplotTop")),
                           tabPanel('Transcript-level analysis', DT::dataTableOutput('TopCisTableTr')),
                           bsModal("ToCisPlotTr", "Top Cis eQTL", "PlotToCisTr", size = "large",plotOutput("eQTLplotTopTr"))
                         )   
                       )
                       
                     )
            ),
           tabPanel("Trans eQTLs",
                    sidebarLayout(
                      sidebarPanel(
                        radioButtons("p_type", "Maximum p-value", c('Uncorrected p-values' = 'nominal_p', 'FDR corrected p-values (q-values)'= 'qvalue'), selected = 'qvalue', inline = FALSE,
                                     width = NULL),
                        sliderInput("pvalue", "p-value:", 
                                    min = 0, max = 1, value = 0.1, step= 0.01),
                        textInput("typedPval", "Type p-value", value=.1),
                        conditionalPanel(
                          'input.trans_dataset === "Gene-level trans analysis"',
                          HTML("<strong>Select row to plot data</strong><br>"),
                          actionButton("PlotTrans", "Plot")
                        ),
                        conditionalPanel(
                          'input.trans_dataset === "Transcript-level trans analysis"',
                          HTML("<strong>Select row to plot data</strong><br>"),
                          actionButton("PlotTransTr", "Plot")
                        )
                      ),
                      mainPanel(
                        tabsetPanel(
                          id = 'trans_dataset',
                          tabPanel('Gene-level trans analysis', DT::dataTableOutput('TransTable')),
                          bsModal("TransPlot", "Trans eQTL", "PlotTrans", size = "large",plotOutput("eQTLplotTrans")),
                          tabPanel('Transcript-level trans analysis', DT::dataTableOutput('TransTableTr')),
                          bsModal("TransPlotTr", "Trans eQTL", "PlotTransTr", size = "large",plotOutput("eQTLplotTransTr"))
                        )   
                      )
                      
                    )
           )
)
