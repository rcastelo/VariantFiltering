library(graph)
library(RBGL)
library(Rgraphviz)
library(shiny)
library(shinyTree)
library(shinyjs)
library(shinythemes)
library(VariantFiltering)
library(DT)

fluidPage(theme = shinytheme("cerulean"),
          shinyjs::useShinyjs(),
          tags$head(
            tags$style(
              HTML(".dataTables_filter, .dataTables_info .shiny-bound-input .shiny-bound-output, .dataTables_length, #page { display: none; }"),
              HTML("#order{margin-top:-20px; margin-bottom:30px;} "),
              HTML("#closeSave{margin-top:20px;} "),
              HTML(".btn-warning { float:right; }")
            )
          ),
          uiOutput("ui")
          
)
