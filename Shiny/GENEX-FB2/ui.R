#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)


# Application title
#titlePanel("Gene Expression in the Fetal Brain: Sex Biases"),

navbarPage("Gene Expression in the Fetal Brain: EQTL:",
            tabPanel("Cis eQTLs",
                     sidebarLayout(
                       sidebarPanel(
                         radioButtons("p_type", "Maximum p-value", c('Uncorrected p-values' = 'pvalue', 'FDR corrected p-values (q-values)'= 'padj'), selected = 'padj', inline = FALSE,
                                      width = NULL),
                         sliderInput("pvalue", textOutput("SciNotation"), 
                                     min = -80, max = -2, value = -10),
                         conditionalPanel(
                           'input.dataset === "Top Cis"',
                           plotOutput("eQTLplotTop", height=200)
                         )
                        ),
                       mainPanel(
                         tabsetPanel(
                           id = 'dataset',
                           tabPanel('Top Cis', DT::dataTableOutput('TopCisTable'))
                         )   
                       )
                       
                     )
            )
)
