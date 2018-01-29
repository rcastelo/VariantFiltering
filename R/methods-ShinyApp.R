browseVariants <- function() {
  ##
  ## Build
  ##
  
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
              shinyjs::disabled(downloadButton("closeSave", strong("Save variants file"), class = "btn-warning"))))
  }
  
  ## validates construction params
  dataValidation <- function(input,output,session) {
    flag <- FALSE
    
    #vcf
    if (!is.null(input$vcfFilename)) {
      vcf <- input$vcfFilename
      flag <- TRUE
    } else flag <- FALSE
    #ped
    if (!is.null(input$pedFilename)) {
      ped <- input$pedFilename
      flag <- TRUE
    } else flag <- FALSE
    
    #orgdb
    if (!is.null(input$orgdb)) {
      orgdb <- input$orgdb
    } else orgdb <- character(0)
    
    #txdb
    if (!is.null(input$txdb)) {
      txdb <- input$txdb
    } else txdb <- character(0)
    
    #snpdb
    if (!is.null(input$snpdb)) {
      snpdb <- input$snpdb
    } else flag <- FALSE
    
    if (!is.null(input$snpdb2)) {
      snpdb2 <- input$snpdb2
    } else flag <- FALSE
    
    #otherAnnotations
    if (!is.null(input$otherAnnotations)) {
      otherAnnotations <- input$otherAnnotations
    } else otherAnnotations <- character(0)
    
    if(flag) {
      sn <- c(snpdb,snpdb2)
      withProgress(message = 'Building', value = 0, {
        incProgress(1, detail = "This may take a while...")
        result = tryCatch({
          vfpardev <<- VariantFilteringParam(vcfFilename = vcf$datapath,pedFilename=ped$datapath,
                                             orgdb=orgdb, txdb=txdb, snpdb=sn, otherAnnotations=otherAnnotations)
        }, warning = function(w) {
          #warning-handler-code
        }, error = function(e) {
          shinyjs::alert("Select correct file!")
        })
      })
    }
  }
  
  ## renderizes construction progress
  renderProgress <- function(input, output, session, progress, message) {
    switch(progress,
           "1" = {
             output$progress <- renderUI(HTML("<div class='progress'><div class='progress-bar' role='progressbar' aria-valuenow='40' 
                                              aria-valuemin='0' aria-valuemax='100' style='width:33%'></div></div>"))
           },
           "2" = {
             output$progress <- renderUI(HTML("<div class='progress'><div class='progress-bar' role='progressbar' aria-valuenow='40' 
                                              aria-valuemin='0' aria-valuemax='100' style='width:66%'></div></div>")) 
           },
           "3" = {
             output$progress <- renderUI(HTML("<div class='progress'><div class='progress-bar' role='progressbar' aria-valuenow='40' 
                                              aria-valuemin='0' aria-valuemax='100' style='width:100%'></div></div>"))
           })
    output$info <- renderText(message)
           }
  
  ##
  ## Browse
  ##
  
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
            if(flt == "SOterms")
            {
              shinyjs::enable("tree")
              active(filters(vfResultsObj))[[flt]] <<- TRUE
              datatable_generation(input,output,session,global)
            }
            else
            {
              shinyjs::enable(flt)
              shinyjs::enable(paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""))
              active(filters(vfResultsObj))[[flt]] <<- TRUE
              datatable_generation(input,output,session,global)
            }
          }
          else
          {
            if(flt == "SOterms")
            {
              shinyjs::disable("tree")
              active(filters(vfResultsObj))[[flt]] <<- FALSE
              datatable_generation(input,output,session,global)
            }
            else
            {
              shinyjs::disable(flt)
              shinyjs::disable(paste(flt,names(cutoffs(vfResultsObj)[[flt]][i]),sep=""))
              active(filters(vfResultsObj))[[flt]] <<- FALSE
              datatable_generation(input,output,session,global)
            }
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
               if(flt == "SOterms") {
                 observeEvent(input$tree,{
                   if(length(unlist(shinyTree::get_selected(input$tree)) == 0)) {
                     change(cutoffs(vfResultsObj), flt) <<- unlist(shinyTree::get_selected(input$tree))                     
                   }
                   else {
                     change(cutoffs(vfResultsObj), flt) <<- "Any"
                   }
                   datatable_generation(input,output,session,global)
                 })
               }
               else {
                 observeEvent(eval(parse(text = paste('input$',flt,sep=''))), {
                   change(cutoffs(vfResultsObj), flt) <<- eval(parse(text = paste("input$",flt,sep="")))
                   datatable_generation(input,output,session,global)
                 })
               }
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
  
  build <- function(input,output,session) {
    #Button comprovation
    if(!exists("vfpardev")) {
      shinyjs::disable("run")
      renderProgress(input,output,session,1, "Step 1: Build your vfParam object by selecting the parameters in the left.")
    } else renderProgress(input,output,session, 2, "Step 2: Run unrelated individuals to build your final object.")
    
    if(!exists("vfResultsObj")) {
      shinyjs::disable("browse")
      shinyjs::disable("closeSave")
    } else renderProgress(input,output,session, 3, "Step 3: Browse the results in our visualitzation app or save the generated object
                       in your workspace.")
    
    #Buttons events
    observeEvent(input$returnButton, {
      dataValidation(input,output,session)
      if(exists("vfpardev")) {
        shinyjs::enable("run")
        renderProgress(input,output,session, 2, "Step 2: Run unrelated individuals to build your final object.")
      }
    })
    
    observeEvent(input$run, {
      if(exists("vfpardev")) {
        withProgress(message = 'Building unrelated individuals', value = 0, {
          incProgress(1, detail = "This may take a while...")
          vfResultsObj <<- unrelatedIndividuals(vfpardev) 
        })
      }
      if(exists("vfResultsObj")) {
        shinyjs::enable("browse")
        shinyjs::enable("closeSave")
        renderProgress(input,output,session, 3, "Step 3: Browse the results in our visualitzation app or save the generated object
                       in your workspace.")
      }
      })
    
    observeEvent(input$browse, {
      remove(vfpardev)
      shinyjs::runjs("location.reload();")
    })
    
    output$closeSave <- downloadHandler(
      filename = function() {
        paste("my_data",".rds", sep="")
      },
      content = function(file){
        saveRDS(vfResultsObj, file)
      })
    }
  
  browse <- function(input,output,session) {
    
    global <- reactiveValues()
    
    generateTree <- function(n,g) {
      sapply(n,
             function(x){
               if(length(unlist(edges(g)[[x]])) == 0) {
                 return(1)
               }else{
                 generateTree(edges(g)[[x]], g)
               }
             }, USE.NAMES = TRUE, simplify = FALSE)
    }
    
    datatable_generation(input,output,session,global)
    
    gn <- sog(vfResultsObj,reverse=TRUE)
    to <- tsort(gn)
    
    tree1 <- generateTree(to[1],gn)
    nodeD <- nodeData(gn, nodes(gn), "label")
    
    ## renames the tree
    rename_items <- function(item){
      if (is.list(item)){
        item <- lapply(item,rename_items)
        names(item) <- unname(nodeD[names(item)])
      }
      item
    }
    
    tree2<-rename_items(tree1)
    
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
    
    observeEvent(input$close, {
      stopApp(vfResultsObj)
    })
    
  }
  
  app <- list(ui = NULL, server = NULL)
  
  app$ui <- fluidPage(theme = shinytheme("cerulean"),
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
  
  app$server <- function(input, output, session) {
    if(!exists("vfResultsObj")) {
      output$ui <- renderUI({
        fluidPage(theme = shinytheme("cerulean"), shinyjs::useShinyjs(),
                  tags$head(
                    tags$style(
                      HTML(".btn-warning { float:right; }")
                    )),
                  fluidRow(
                    conditionalPanel(id= "buildPanel",
                                     column(6),
                                     mainDataInput("mainInput"),
                                     column(3,
                                            h3("Progress"),
                                            htmlOutput("progress"),
                                            textOutput("info"))
                    )
                  )
        )
      })
      build(input,output,session)  
    }
    else {
      output$ui <- renderUI({
        inputData("inputData", vfResultsObj) 
      })
      browse(input,output,session)
    }
  }
  runApp(app)
}