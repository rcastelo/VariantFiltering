library(graph)
library(RBGL)
library(Rgraphviz)
library(shiny)
library(shinyTree)
library(shinyjs)
library(shinythemes)
library(VariantFiltering)
library(DT)

## build construction params
mainDataInput <- function(id) {
  ns <- NS(id)
  mainPanel(width = 6,
            h2("Select data source:"),
            wellPanel(fluidRow(
              column(
                12,
                #vcfFilename
                fileInput("vcfFilename", strong("Choose vcf:"),
                          accept = c(".vcf",".vcf.gz",".vcf.bgz")
                ),
                #pedFilename
                fileInput("pedFilename", strong("Choose ped:"),
                          accept = c(".ped")
                ),
                HTML("<br>"),
                #orgdb
                selectInput("orgdb", strong("Select orgdb:"),
                            VariantFiltering:::.orgdbpkgs()),
                HTML("<br>"),
                #txdb
                selectInput("txdb", strong("Select txdb:"),
                            VariantFiltering:::.txdbpkgs()),
                HTML("<br>"),
                #snpdb
                selectInput("snpdb", strong("Select snpdb:"),
                            VariantFiltering:::.snpdbpkgs()),
                selectInput("snpdb2","",
                            VariantFiltering:::.xtrasnpdbpkgs()),
                HTML("<br>"),
                #otherAnnotation
                checkboxGroupInput("otherAnnotations",strong("Select other annotations:"),VariantFiltering:::.otherannotations())
                
              )
            ),
            HTML("<br>"),
            actionButton("returnButton", strong("Build"), class = "btn-primary"),
            shinyjs::disabled(actionButton("run", strong("Run"), class = "btn-info")),
            shinyjs::disabled(actionButton("browse", strong("Browse Variants"), class = "btn-success")),
            shinyjs::disabled(downloadButton("closeSave", strong("Save parameters file"), class = "btn-warning"))))
}

## filtering control of logical values
logicalPanel <- function(flt,tab) {
  conditionalPanel(condition=paste("input.tsp ==", tab, sep = ""), 
                   checkboxInput(paste("active",flt,sep=""), strong(flt),
                                 active(filters(vfResultsObj))[[flt]]),
                   hidden(p(id=paste0("info",flt),tags$i(filtersMetadata(vfResultsObj)[!is.na(filtersMetadata(vfResultsObj)$AnnoGroup),"Description"][[flt]]))),
                   checkboxGroupInput(flt,flt,c(names(cutoffs(vfResultsObj)[[flt]])), selected = names(which(c(cutoffs(vfResultsObj)[[flt]])))))
}

## filtering control of a vector of character strings
characterVectorPanel <- function(flt,tab) {
  conditionalPanel(condition=paste("input.tsp ==", tab, sep = ""),
                   checkboxInput(paste("active",flt,sep=""), strong(flt),
                                 active(filters(vfResultsObj))[[flt]]),
                   hidden(p(id=paste0("info",flt),tags$i(filtersMetadata(vfResultsObj)[!is.na(filtersMetadata(vfResultsObj)$AnnoGroup),"Description"][[flt]]))),
                   selectInput(flt, flt,unname(cutoffs(vfResultsObj)[[flt]])))
}

## filtering control of a character string 
characterPanel <- function(flt,tab) {
  conditionalPanel(condition=paste("input.tsp ==", tab, sep = ""),
                   checkboxInput(paste("active",flt,sep=""), strong(flt),
                                 active(filters(vfResultsObj))[[flt]]),
                   hidden(p(id=paste0("info",flt),tags$i(filtersMetadata(vfResultsObj)[!is.na(filtersMetadata(vfResultsObj)$AnnoGroup),"Description"][[flt]]))),
                   #textInput(flt, flt, value = cutoffs(vfResultsObj)[[flt]][1])
                   shinyTree::shinyTree("tree", flt, checkbox = TRUE))
}

## filtering control of a vector of numeric values
numericVectorPanel <- function(flt,tab) {
  conditionalPanel(condition=paste("input.tsp ==", tab, sep = ""),
                   checkboxInput(paste("active",flt,sep=""), strong(flt),
                                 active(filters(vfResultsObj))[[flt]]),
                   hidden(p(id=paste0("info",flt),tags$i(filtersMetadata(vfResultsObj)[!is.na(filtersMetadata(vfResultsObj)$AnnoGroup),"Description"][[flt]]))),
                   selectInput(flt, flt,
                               c(cutoffs(vfResultsObj)[[flt]])))
}

## filtering control of a numeric value
numericPanel <- function(flt,tab) {
  conditionalPanel(condition=paste("input.tsp ==", tab, sep = ""),
                   checkboxInput(paste("active",flt,sep=""), strong(flt),
                                 active(filters(vfResultsObj))[[flt]]),
                   hidden(p(id=paste0("info",flt),tags$i(filtersMetadata(vfResultsObj)[!is.na(filtersMetadata(vfResultsObj)$AnnoGroup),"Description"][[flt]]))),
                   numericInput(flt, flt, value = cutoffs(vfResultsObj)[[flt]][1], min = 0, max = 0.5))
}

## filtering control for a factor value
factorPanel <- function(flt,tab) {
  conditionalPanel(condition=paste("input.tsp ==", tab, sep = ""),
                   checkboxInput(paste("active",flt,sep=""), strong(flt),
                                 active(filters(vfResultsObj))[[flt]]),
                   hidden(p(id=paste0("info",flt),tags$i(filtersMetadata(vfResultsObj)[!is.na(filtersMetadata(vfResultsObj)$AnnoGroup),"Description"][[flt]]))),
                   selectInput(flt, flt,levels(cutoffs(vfResultsObj)[[flt]]), selected = as.character(cutoffs(vfResultsObj)[[flt]])))
}

## filtering control of a list of values of different classes

listPanel <- function(flt,tab) {
  conditionalPanel(condition=paste("input.tsp ==", tab, sep = ""), checkboxInput(paste("active",flt,sep=""), strong(flt),
                                                                                 active(filters(vfResultsObj))[[flt]]),
                   hidden(p(id=paste0("info",flt),tags$i(filtersMetadata(vfResultsObj)[!is.na(filtersMetadata(vfResultsObj)$AnnoGroup),"Description"][[flt]]))),
                   lapply(1:NROW(cutoffs(vfResultsObj)[[flt]]), function(i)
                     if(class(cutoffs(vfResultsObj)[[flt]][[i]]) == "numeric")
                     {
                       if(length(cutoffs(vfResultsObj)[[flt]][[i]]) > 1)
                       {
                         selectInput(paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""), names(cutoffs(vfResultsObj)[[flt]][i]),
                                     c(cutoffs(vfResultsObj)[[flt]]))
                       }
                       else
                       {
                         numericInput(paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""), names(cutoffs(vfResultsObj)[[flt]][i]), 
                                      value = cutoffs(vfResultsObj)[[flt]][[i]], min = 0, max = 0.5)
                       }
                     }
                     else if(class(cutoffs(vfResultsObj)[[flt]][[i]]) == "logical")
                     {
                       checkboxGroupInput(paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""), names(cutoffs(vfResultsObj)[["maxMAF"]][i]),
                                          c(names(cutoffs(vfResultsObj)[[flt]][[i]])),selected = names(which(c(cutoffs(vfResultsObj)[[flt]][[i]]))))
                     }
                     else if(class(cutoffs(vfResultsObj)[[flt]][[i]]) == "character")
                     {
                       if(length(cutoffs(vfResultsObj)[[flt]][[i]]) > 1)
                       {
                         selectInput(paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""), names(cutoffs(vfResultsObj)[[flt]][i]),
                                     c(cutoffs(vfResultsObj)[[flt]][[i]]))
                       }
                       else
                       {
                         ##Selected
                         textInput(paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""), names(cutoffs(vfResultsObj)[[flt]][i]), 
                                   value = cutoffs(vfResultsObj)[[flt]][i])
                       }
                     }
                     else if(class(cutoffs(vfResultsObj)[[flt]][[i]]) == "factor")
                     {
                       selectInput(paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""), names(cutoffs(vfResultsObj)[[flt]][i]),
                                   levels(cutoffs(vfResultsObj)[[flt]][[i]]),selected = as.character(cutoffs(vfResultsObj)[[flt]][[i]]))
                     }
                     else
                     {
                       ##Variant list without cutoffs
                     }
                   )
  )
}

## build the filtering control panel
generateFiltering <- function(vfResultsObj) {
  lapply(1:NROW(filtersMetadata(vfResultsObj)[!is.na(filtersMetadata(vfResultsObj)$AnnoGroup),]), function(j) {
    tab <- tolower(paste("'","'", sep=filtersMetadata(vfResultsObj)[!is.na(filtersMetadata(vfResultsObj)$AnnoGroup),][,2][j]))
    flt <- rownames(filtersMetadata(vfResultsObj)[!is.na(filtersMetadata(vfResultsObj)$AnnoGroup),])[j]
    switch(class(cutoffs(vfResultsObj)[[flt]]),
           "logical" = {
             logicalPanel(flt,tab)
           },
           "character" = {
             if(length(cutoffs(vfResultsObj)[[flt]]) > 1) {
               characterVectorPanel(flt,tab)
             }
             else {
               characterPanel(flt,tab)
             }
           },
           "numeric" = {
             if(length(cutoffs(vfResultsObj)[[flt]]) > 1) {
               numericVectorPanel(flt,tab)
             }
             else {
               numericPanel(flt,tab)
             }
           },
           "CutoffsList" = {
             listPanel(flt,tab)
           },
           "factor" = {
             factorPanel(flt,tab)
           },
           NULL = {
             conditionalPanel(condition=paste("input.tsp ==", tab, sep = ""), checkboxInput(paste("active",flt,sep=""), 
                                                                                            strong(flt),active(filters(vfResultsObj))[[flt]]),
                              hidden(p(id=paste0("info",flt),tags$i(filtersMetadata(vfResultsObj)[!is.na(filtersMetadata(vfResultsObj)$AnnoGroup),"Description"][[flt]]))))
           }
    )
  })
}

## layout the UI of the controls
inputData <- function(id, vfResultsObj) {
  ns <- NS(id)
  
  tabs <- lapply(names(annoGroups(vfResultsObj)),
                 function(tabname) tabPanel(tabname,
                                            DT::dataTableOutput(paste0("table", tabname)),
                                            value=tolower(tabname)))
  tabs[[length(tabs)+1]] <- "tsp"
  names(tabs) <- c(rep("", length(tabs)-1), "id")
  
  fluidRow(
    #UI declaration
    column(
      3,
      h2("Select data to filter:"),
      #Data
      wellPanel(fluidRow(
        column(
          12,
          selectInput("sort", "Select column to sort by:",levels(sortings(vfResultsObj)[["criterion"]]), selected = as.character(sortings(vfResultsObj)[["criterion"]])),
          radioButtons("order","",choiceNames = list("Decreasing", "Increasing"), choiceValues = list(FALSE, TRUE), selected = sortings(vfResultsObj)[["decreasing"]]),
          generateFiltering(vfResultsObj),
          actionButton("close", "Save and close", class = "btn btn-warning")
        )
      ))
    ),
    column(1, uiOutput("pageFilter")),
    column(8, mainPanel(h2("Variant Visualization"),do.call(tabsetPanel, tabs)))
  )
}

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