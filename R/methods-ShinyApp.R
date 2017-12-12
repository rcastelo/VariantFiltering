browseVariants <- function(vfResultsObj) {

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
                   textInput(flt, flt, value = cutoffs(vfResultsObj)[[flt]][1]))
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
inputData <- function(id, vfrobj) {
  ns <- NS(id)

  tabs <- lapply(names(annoGroups(vfrobj)),
                 function(tabname) tabPanel(tabname,
                                            DT::dataTableOutput(paste0("table", tabname)),
                                            value=tolower(tabname)))
  tabs[[length(tabs)+1]] <- "tsp"
  names(tabs) <- c(rep("", length(tabs)-1), "id")

  fluidRow(
    #UI declaration
    column(
      3,
      #Data
      wellPanel(fluidRow(
        column(
          12,
          selectInput("sort", "Select column to sort by:",levels(sortings(vfResultsObj)[["criterion"]]), selected = as.character(sortings(vfResultsObj)[["criterion"]])),
          radioButtons("order","",choiceNames = list("Decreasing", "Increasing"), choiceValues = list(FALSE, TRUE), selected = sortings(vfResultsObj)[["decreasing"]]),
          generateFiltering(vfResultsObj),
          actionButton("closeSave", "Save and Close")
        )
      ))
    ),
    column(1, uiOutput("pageFilter")),
    column(8, mainPanel(do.call(tabsetPanel, tabs)))
  )
}

## populate the table of variants
datatable_generation <- function(input, output, session, global) {
  
  lapply(names(annoGroups(vfResultsObj)), function(tabname){
    observe({global$val <- eval(parse(text = paste0("input$table",tabname,"_state$start"))) / eval(parse(text = paste0("input$table",tabname,"_state$length"))) + 1})
  })
  
  output$pageFilter <- renderUI({
    numericInput("page", "", global$val, min = 1)
  })
  
  filteredVariantsReact <- reactive({
    vdf <- data.frame()
    fvxsam <- filteredVariants(vfResultsObj)
    fv <- fvxsam[[1]]
    vdf <- DataFrame(VarID=fv$VARID, CHR=seqnames(fv),
                     POS=start(fv), POSITION=start(fv),
                     mcols(fv))
    vdf$REF <- ref(fv)
    vdf$ALT <- alt(fv)
    vdf <- as.data.frame(vdf)
    rownames(vdf) <- NULL
    vdf
  })
  
  withProgress(message = "Loading...",
               lapply(names(annoGroups(vfResultsObj)), function(tabname) {
                 output[[paste0("table",tabname)]] <- DT::renderDataTable({
                   DT::datatable(filteredVariantsReact()[,annoGroups(vfResultsObj)[[tabname]]], options = list(ordering=F, pageLength = 10, stateSave = TRUE))
                 })
               })
  )
  
}

## activate filters
filtersActivation <- function(input, output, session, global) {
  lapply(rownames(filtersMetadata(vfResultsObj)[!is.na(filtersMetadata(vfResultsObj)$AnnoGroup),]), function(flt)
    lapply(1:NROW(cutoffs(vfResultsObj)[[flt]]), function(i)
      observeEvent(eval(parse(text = paste("input$active",flt,sep=""))), {
        if(eval(parse(text = paste("input$active",flt,sep=""))) == TRUE)
        {
          shinyjs::enable(flt)
          shinyjs::enable(paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""))
          active(filters(vfResultsObj))[[flt]] <<- TRUE
          datatable_generation(input,output,session,global)
        }
        else
        {
          shinyjs::disable(flt)
          shinyjs::disable(paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""))
          active(filters(vfResultsObj))[[flt]] <<- FALSE
          datatable_generation(input,output,session,global)
        }
      })
    )
  )
}

## modify cutoffs

cutoffsModification <- function(input, output, session, global) {
  lapply(rownames(filtersMetadata(vfResultsObj)[!is.na(filtersMetadata(vfResultsObj)$AnnoGroup),]), function(flt)
    switch(class(cutoffs(vfResultsObj)[[flt]]),
     "factor" = {
       observeEvent(eval(parse(text = paste("input$",flt,sep=""))), {
         change(cutoffs(vfResultsObj), flt) <<- eval(parse(text = paste("input$",flt,sep="")))
         datatable_generation(input,output,session,global)
       })
     },
     "numeric" = {
       observeEvent(eval(parse(text = paste('input$',flt,sep=''))), {
         change(cutoffs(vfResultsObj), flt) <<- as.numeric(eval(parse(text = paste("input$",flt,sep=""))))
         datatable_generation(input,output,session,global)
       })
     },
     "character" = {
       observeEvent(eval(parse(text = paste('input$',flt,sep=''))), {
         change(cutoffs(vfResultsObj), flt) <<- eval(parse(text = paste("input$",flt,sep="")))
         datatable_generation(input,output,session,global)
       })
     },
     "logical" = {
       observeEvent(eval(parse(text = paste('input$',flt,sep=''))), {
         change(cutoffs(vfResultsObj),flt) <<- setNames(names(cutoffs(vfResultsObj)[[flt]]) %in% eval(parse(text = paste('input$',flt,sep=''))),
                                                        names(cutoffs(vfResultsObj)[[flt]]))
         datatable_generation(input,output,session,global)
       })
     },
     "CutoffsList" = {
       lapply(1:NROW(cutoffs(vfResultsObj)[[flt]]), function(i) {
         switch(class(cutoffs(vfResultsObj)[[flt]][[i]]),
          "numeric" = {
            observeEvent(eval(parse(text = paste('input$',paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""),sep=''))), {
              change(cutoffs(vfResultsObj)[[flt]], names(cutoffs(vfResultsObj)[[flt]][i])) <<- as.numeric(eval(parse(text = paste('input$',paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""),sep=''))))
              datatable_generation(input,output,session,global)
            })
          },
          "logical" = {
            observeEvent(eval(parse(text = paste('input$',paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""),sep=''))), {
              change(cutoffs(vfResultsObj)[[flt]], names(cutoffs(vfResultsObj)[[flt]][i])) <<- setNames(names(cutoffs(vfResultsObj)[[flt]][[i]]) %in% eval(parse(text = paste('input$',paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""),sep=''))),
                                                                                                        names(cutoffs(vfResultsObj)[[flt]][[i]]))
              datatable_generation(input,output,session,global)
            })
          },
          "character" = {
            observeEvent(eval(parse(text = paste('input$',paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""),sep=''))), {
              change(cutoffs(vfResultsObj)[[flt]], names(cutoffs(vfResultsObj)[[flt]][i])) <<- eval(parse(text = paste('input$',paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""),sep='')))
              datatable_generation(input,output,session,global)
            })
          },
          "factor" = {
            observeEvent(eval(parse(text = paste('input$',paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""),sep=''))), {
              change(cutoffs(vfResultsObj)[[flt]], names(cutoffs(vfResultsObj)[[flt]][i])) <<- eval(parse(text = paste('input$',paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""),sep='')))
              datatable_generation(input,output,session,global)
            })
          }
         )
       })
     }
    )
  )
}

## build the shiny app
app <- list(ui = NULL, server = NULL)

app$ui <- fluidPage(theme = shinytheme("simplex"),
  shinyjs::useShinyjs(),
  tags$head(
    tags$style(
      HTML(".dataTables_filter, .dataTables_info .shiny-bound-input .shiny-bound-output, .dataTables_length, #page { display: none; }"),
      HTML("#order{margin-top:-20px; margin-bottom:30px;} "),
      HTML("#closeSave{margin-top:20px;} ")
    )
  ),
  inputData("inputData", vfResultsObj)
)

app$server <- function(input, output, session) {

  global <- reactiveValues()

  datatable_generation(input,output,session,global)

  #Page filter
  observeEvent({input$page; input$tsp}, {
    lapply(names(annoGroups(vfResultsObj)), function(tabname)
      selectPage(dataTableProxy(paste0("table", tabname)), global$val)
      ## dataTableProxy(paste0("table", tabname)) %>% selectPage(global$val)
    )
  })

  #Activating and deactivating filters
  filtersActivation(input,output,session,global)

  #Changing cutoff values
  cutoffsModification(input,output,session,global)

  #Showing Variant Description
  lapply(rownames(filtersMetadata(vfResultsObj)[!is.na(filtersMetadata(vfResultsObj)$AnnoGroup),]), function(flt) {
    shinyjs::onevent("mouseenter", paste0("active",flt), shinyjs::show(paste0("info",flt)))
    onevent("mouseleave", paste0("active",flt), shinyjs::hide(paste0("info",flt)))
  })

  #Sorting
  observeEvent(input$sort, {
    change(sortings(vfResultsObj),"criterion") <<- input$sort
    datatable_generation(input,output,session,global)
  })
  observeEvent(input$order, {
    change(sortings(vfResultsObj),"decreasing") <<- eval(parse(text = input$order))
    datatable_generation(input,output,session,global)
  })

  #Observe the Save & Close button
  observeEvent(input$closeSave, {
    stopApp(vfResultsObj) #Stops the app and returns the generated_gsva object
  })
}

runApp(app)
}
