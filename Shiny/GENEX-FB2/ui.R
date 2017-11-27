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

navbarPage("Gene Expression in the Fetal Brain: EQTL:",
            tabPanel("Cis eQTLs",
                     sidebarLayout(
                       sidebarPanel(
                         radioButtons("p_type", "Maximum p-value", c('Uncorrected p-values' = 'nominal_p', 'FDR corrected p-values (q-values)'= 'qvalue'), selected = 'qvalue', inline = FALSE,
                                      width = NULL),
                         sliderInput("pvalue", "p-value:", 
                                     min = 0, max = 1, value = 0.1, step= 0.01),
                         textInput("typedPval", "Type p-value", value=.1),
                         conditionalPanel(
                           'input.dataset === "Top Cis"',
                           HTML("<strong>Select row to plot data</strong><br>"),
                           actionButton("PlotToCis", "Plot")
                         )
                        ),
                       mainPanel(
                         tabsetPanel(
                           id = 'dataset',
                           tabPanel('Top Cis', DT::dataTableOutput('TopCisTable')),
                           bsModal("ToCisPlot", "Top Cis eQTL", "PlotToCis", size = "large",plotOutput("eQTLplotTop"))
                         )   
                       )
                       
                     )
            )
)
