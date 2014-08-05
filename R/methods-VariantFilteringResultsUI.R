setMethod("show", signature(object="VariantFilteringResultsUI"),
          function(object) {
            cat("\nVariantFiltering results object\n\n")
            cat(sprintf("  Variants present in a group of %s \n", inheritanceModel(object)))
            cat("  Functional annotation filters\n")
            if (is.na(dbSNPpresent(object)))
                cat(sprintf("    No filtering on presence in %s %s\n", provider(param(object)$snpdb),
                            releaseName(param(object)$snpdb)))
            else
                cat(sprintf("    Present in %s %s: %s\n", provider(param(object)$snpdb),
                            releaseName(param(object)$snpdb), dbSNPpresent(object)))
            cat(sprintf("    Variant type: %s\n", variantType(object)))
            cat(sprintf("    Amino acid change type: %s\n", aaChangeType(object)))
            cat(sprintf("    Populations used for MAF filtering: %s\n", paste(names(MAFpop(object))[MAFpop(object)], collapse=", ")))
            cat(sprintf("    Include MAF NA values: %s\n", ifelse(naMAF(object), "yes", "no")))
            cat(sprintf("    Maximum MAF: %.2f\n", maxMAF(object)))
            whPhCdb <- match("PhastConsDb", sapply(param(object)$otherAnnotations, class))
            if (!is.na(whPhCdb)) {
              if (is.na(minPhastCons(object)))
                cat("    No filtering on nucleotide conservation\n")
              else
                cat(sprintf("    Minimum score for phastCons nucleotide conservation: %.2f\n", minPhastCons(object)))
            }
            whGPSdb <- match("GenePhylostrataDb", sapply(param(object)$otherAnnotations, class))
            if (!is.na(whGPSdb)) {
              if (is.na(minPhylostratum(object)))
                cat("    No filtering on gene conservation\n")
              else
                cat(sprintf("    Minimum conserved gene phylostratum: %s (%d/%d)\n",
                            genePhylostrata(param(object)$otherAnnotations[[whGPSdb]])$Description[minPhylostratum(object)],
                            minPhylostratum(object),
                            nrow(genePhylostrata(param(object)$otherAnnotations[[whGPSdb]]))))
            }
            if (is.na(minCRYP5ss(object)))
              cat("    No filtering on cryptic 5'ss\n")
            else
              cat(sprintf("    Minimum score for cryptic 5'ss: %.2f\n", minCRYP5ss(object)))
            if (is.na(minCRYP3ss(object)))
              cat("    No filtering on cryptic 3'ss\n")
            else
              cat(sprintf("    Minimum score for cryptic 3'ss: %.2f\n", minCRYP3ss(object)))
            vdf <- filteredVariants(object)
            total <- length(vdf)
            vxloc <- table(vdf$LOCATION)
            vxcon <- table(vdf$CONSEQUENCE)
            cat(sprintf("\n  Total number of variants: %d\n", total))
            paddgts <- sprintf("%%%dd", nchar(as.character(total)))
            if (!is.na(match("nonsynonymous", names(vxcon))))
              cat(sprintf("    %s (%.1f%%) are coding non-synonymous\n", sprintf(paddgts, vxcon["nonsynonymous"]), 100*vxcon["nonsynonymous"]/total))
            if (!is.na(match("synonymous", names(vxcon))))
              cat(sprintf("    %s (%.1f%%) are coding synonymous\n", sprintf(paddgts, vxcon["synonymous"]), 100*vxcon["synonymous"]/total))
            if (!is.na(match("nonsense", names(vxcon))))
              cat(sprintf("    %s (%.1f%%) are coding nonsense\n", sprintf(paddgts, vxcon["nonsense"]), 100*vxcon["nonsense"]/total))
            if (!is.na(match("spliceSite", names(vxloc))))
              cat(sprintf("    %s (%.1f%%) located in known splice sites\n", sprintf(paddgts, vxloc["spliceSite"]), 100*vxloc["spliceSite"]/total))
            if (!is.na(match("promoter", names(vxloc))))
              cat(sprintf("    %s (%.1f%%) located in promoter regions\n", sprintf(paddgts, vxloc["promoter"]), 100*vxloc["promoter"]/total))
            if (!is.na(match("fiveUTR", names(vxloc))))
              cat(sprintf("    %s (%.1f%%) located in 5' UTR regions\n", sprintf(paddgts, vxloc["fiveUTR"]), 100*vxloc["fiveUTR"]/total))
            if (!is.na(match("threeUTR", names(vxloc))))
              cat(sprintf("    %s (%.1f%%) located in 3' UTR regions\n", sprintf(paddgts, vxloc["threeUTR"]), 100*vxloc["threeUTR"]/total))
            if (!is.na(match("intron", names(vxloc))))
              cat(sprintf("    %s (%.1f%%) located in intronic regions\n", sprintf(paddgts, vxloc["intron"]), 100*vxloc["intron"]/total))
          })


setMethod("selectIndividual", signature(x="VariantFilteringResultsUI"),
          function(x) {
              x@indselected
          })

setReplaceMethod("selectIndividual", signature(x="VariantFilteringResultsUI", value="ANY"),
                 function(x, value) {
                     x@indselected <- as.character(value)
                     x
          })

setMethod("selectGene", signature(x="VariantFilteringResultsUI"),
          function(x) {
              x@selectgene
          })

setReplaceMethod("selectGene", signature(x="VariantFilteringResultsUI", value="ANY"),
                 function(x, value) {
                     x@selectgene <- as.character(value)
                     x
          })


## get variants after applying all filters
setMethod("filteredVariants", signature(x="VariantFilteringResultsUI"), 
          function(x, unusedColumns.rm=FALSE) {
            vars <- allVariants(x)
            rowsMask <- rep(TRUE, length(vars))

            ## selection of individuals
            if (length(selectIndividual(x)) > 1) {
              if (any(is.na(selectIndividual(x)))) {
                stop("No NAs allowed in the vector of individuals.")
              }

              hom <- sapply(selectIndividual(x), function(z) sum(which(vars$HomozygousAlt == z)) > 0)
              maskIndHom <- sapply(seq(nrow(hom)), function(z) any(hom[z, ] == TRUE))

              het <- sapply(selectIndividual(x), function(z) sum(which(vars$Heterozygous == z)) > 0)
              maskIndHet <- sapply(seq(nrow(het)), function(z) any(het[z, ] == TRUE))

              rowsMask <- rowsMask & (maskIndHet | maskIndHom)

            } else {

              if (!is.na(selectIndividual(x))) {
                maskIndHom <- sum(which(vars$HomozygousAlt == selectIndividual(x))) > 0
                maskIndHet <- sum(which(vars$Heterozygous == selectIndividual(x))) > 0
                rowsMask <- rowsMask & (maskIndHet | maskIndHom)
              }

            }

            ## select concrete gene
            if (!is.na(selectGene(x))) {
                maskGene <- vars$GENE == selectGene(x)
                maskGene[is.na(maskGene)] <- FALSE
                rowsMask <- rowsMask & maskGene    
            }
            
            ## presence in dbSNP
            if (!is.na(dbSNPpresent(x))) {
              maskNAdbSNP <- is.na(vars$dbSNP)
              if (dbSNPpresent(x) == "Yes")
                rowsMask <- rowsMask & !maskNAdbSNP
              else
                rowsMask <- rowsMask & maskNAdbSNP
            }

            ## presence in OMIM
            if (!is.na(OMIMpresent(x))) {
              maskNAomim <- is.na(vars$OMIM)
              if (OMIMpresent(x) == "Yes")
                rowsMask <- rowsMask & !maskNAomim
              else
                rowsMask <- rowsMask & maskNAomim
            }

            ## type of variant
            if (variantType(x) != "Any")
              rowsMask <- rowsMask & vars$TYPE == variantType(x)

            ## type of amino acid change
            if (aaChangeType(x) != "Any")
              rowsMask <- rowsMask & vars$AAchangeType == aaChangeType(x)

            ## minimum allele frequency
            mtNoMAF <- NULL
            if (!is.na(match("MafDb", sapply(param(x)$otherAnnotations, class)))) {
              vars$maxMAF <- rep(NA_real_, length(vars))
              if (any(MAFpop(x))) {
                vars$maxMAF <- do.call(pmax, c(as.list(mcols(vars[, names(MAFpop(x))[MAFpop(x)]])), na.rm=TRUE))
                naMAFmask <- rep(TRUE, length(vars))
                if (naMAF(x))
                  vars$maxMAF[is.na(vars$maxMAF)] <- -Inf
                else
                  vars$maxMAF[is.na(vars$maxMAF)] <- Inf

                rowsMask <- rowsMask & naMAFmask & vars$maxMAF <= maxMAF(x)
                rowsMask[is.na(rowsMask)] <- FALSE

                vars$maxMAF[!is.finite(vars$maxMAF)] <- NA_real_
                if (any(!MAFpop(x)))
                  mtNoMAF <- match(names(MAFpop(x))[!MAFpop(x)], colnames(mcols(vars)))
              } else {
                mtNoMAF <- match(names(MAFpop(x)), colnames(mcols(vars)))
                if (!naMAF(x))
                  rowsMask <- rep(FALSE, length(vars))
              }
            }

            ## nucleotide conservation
            mtNoMinPhastCons <- NULL
            if (!is.na(match("PhastConsDb", sapply(param(x)$otherAnnotations, class)))) {
              if (is.na(minPhastCons(x)))
                mtNoMinPhastCons <- match("phastCons", colnames(mcols(vars)))
              else {
                rowsMask <- rowsMask & vars$phastCons >= minPhastCons(x)
                rowsMask[is.na(rowsMask)] <- FALSE
              }
            }

            ## gene conservation
            mtNoMinPhylostratum <- NULL
            if (!is.na(match("GenePhylostrataDb", sapply(param(x)$otherAnnotations, class)))) {
              if (is.na(minPhylostratum(x)))
                mtNoMinPhylostratum <- grep("GenePhylostratum", colnames(mcols(vars)))
              else {
                rowsMask <- rowsMask & vars$GenePhylostratumIndex <= minPhylostratum(x)
                rowsMask[is.na(rowsMask)] <- FALSE
              }
            }

            ## if any of the 5' cryptic ss or 3' cryptic ss meet the cutoff, select the row
            crypssMask <- rep(FALSE, length(vars))
            mtNoCRYP5ss <- NULL
            if (is.na(minCRYP5ss(x)))
              mtNoCRYP5ss <- grep("CRYP5ss", colnames(mcols(vars)))
            else {
              cryp5ssMask <- vars$CRYP5ssALT >= minCRYP5ss(x)
              cryp5ssMask[is.na(cryp5ssMask)] <- FALSE
              crypssMask <- crypssMask | cryp5ssMask
            }
            mtNoCRYP3ss <- NULL
            if (is.na(minCRYP3ss(x)))
              mtNoCRYP3ss <- grep("CRYP3ss", colnames(mcols(vars)))
            else {
              cryp3ssMask <- vars$CRYP3ssALT >= minCRYP3ss(x)
              cryp3ssMask[is.na(cryp3ssMask)] <- FALSE
              crypssMask <- crypssMask | cryp3ssMask
            }
            ## if no filter on 5' and 3' cryptic ss is set, then select all rows
            if (is.na(minCRYP5ss(x)) && is.na(minCRYP3ss(x)))
              crypssMask <- rep(TRUE, length(vars))

            rowsMask <- rowsMask & crypssMask

            colsMask <- setdiff(1:ncol(mcols(vars)), mtNoMAF)

            if (unusedColumns.rm) ## remove data columns that are not used for filtering
              colsMask <- setdiff(colsMask, c(mtNoMinPhastCons, mtNoMinPhylostratum, mtNoCRYP5ss, mtNoCRYP3ss))

            vars[rowsMask, colsMask]
          })




             ##########################
             ####### shiny app ########
             ##########################


setMethod("reportVariants", signature(vfResultsObj="VariantFilteringResultsUI"),
          function(vfResultsObj, type=c("shiny", "csv", "tsv"), file=NULL) {
              
  type <- match.arg(type)

  if (class(vfResultsObj) != "VariantFilteringResultsUI")
    stop("Input argument 'vfResultsObj' should be an object of class 'VariantFilteringResultsUI'.")

  if ((type == "csv" || type == "tsv") && is.null(file))
    stop("If type=\"csv\" or type=\"tsv\" then the input argument 'file' cannot be NULL.")

  if (type == "csv" || type == "tsv") {
    varsdf <- filteredVariants(vfResultsObj)
    varsdf <- as.data.frame(DataFrame(VarID=names(varsdf), CHR=seqnames(varsdf),
                                      POS=start(varsdf), mcols(varsdf)))
    if (type == "csv")
      write.csv(varsdf, file=file, quote=FALSE, row.names=FALSE, col.names=TRUE)
    else
      write.table(varsdf, file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
   

    bn <- basename(file)
    message(sprintf("Written %s", bn))
    fullpath <- file
    if (dirname(fullpath) == ".")
      fullpath <- file.path(getwd(), bn)
    return(fullpath)
  }
  
  ######## parameters for the UCSC genome browser, by now we have to hard-code "Hsapiens"
  species <- "Hsapiens"
  genomeVersion <- as.vector(genome(allVariants(vfResultsObj))[1])

  indlist <- gsub("homALT_", "", colnames(mcols(filteredVariants(vfResultsObj)))[grep("homALT_", colnames(mcols(filteredVariants(vfResultsObj))))])
  indlist <- sort(indlist)

  genelist <- unique(allVariants(vfResultsObj)$GENE)

  annotationObjClasses <- sapply(param(vfResultsObj)$otherAnnotations, class)
  mtMafDb <- match("MafDb", annotationObjClasses)
  mtPhastConsDb <- match("PhastConsDb", annotationObjClasses)
  mtGenePhylostrataDb <- match("GenePhylostrataDb", annotationObjClasses)
  phylostrata <- "Unavailable"
  if (!is.na(mtGenePhylostrataDb))
    phylostrata <- rev(genePhylostrata(param(vfResultsObj)$otherAnnotations[[mtGenePhylostrataDb]])$Description)
  
  app <- list(ui=NULL, server=NULL)
  app$ui <- pageWithSidebar(
              headerPanel('R/BioC VariantFiltering Web App'),
              sidebarPanel(
                ## genome tab
                conditionalPanel(condition="input.tsp == 'genome'",
                                 checkboxInput('individualFlag',
                                               'Filter by individual', FALSE)), 
                conditionalPanel(condition="input.tsp == 'genome'&& input.individualFlag == true",
                                 selectInput("selectIndividual", "Individual:",
                                             choices=indlist)),
                conditionalPanel(condition="input.tsp == 'genome'",
                                 checkboxInput('dbSNPpresentFlag',
                                               'Filter by presence in dbSNP', FALSE)),
                conditionalPanel(condition="input.tsp == 'genome' && input.dbSNPpresentFlag == true",
                                 selectInput("dbSNPpresent", "Present in dbSNP:",
                                             choices=c("Yes", "No"))),
                conditionalPanel(condition="input.tsp == 'genome'", selectInput("variantType", "Variant Type:",
                                                                    choices=c("Any", "SNV", "InDel", "MNV"))),
                ## gene tab
                conditionalPanel(condition="input.tsp == 'gene'",
                                 checkboxInput('geneFlag',
                                               'Filter by gene', FALSE)), 
                conditionalPanel(condition="input.tsp == 'gene'&& input.geneFlag == true",
                                 selectInput("selectGene", "Gene:",
                                             choices=genelist)),
                 conditionalPanel(condition="input.tsp == 'gene'",
                                 checkboxInput('OMIMpresentFlag',
                                               'Filter by presence in OMIM', FALSE)),
                conditionalPanel(condition="input.tsp == 'gene' && input.OMIMpresentFlag == true",
                                 selectInput("OMIMpresent", "Present in OMIM:",
                                             choices=c("Yes", "No"))),
                ## protein tab
                conditionalPanel(condition="input.tsp == 'protein'", selectInput("aaChangeType", "Amino acid change type:",
                                                                     choices=c("Any", "Radical", "Conservative"))),
                ## MAF tab
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('naMAF', 'Keep variants without MAF', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", numericInput('maxMAF', 'Maximum MAF:', 1.00)),
                conditionalPanel(condition="input.tsp == 'maf'", helpText("Note: the maximum MAF cutoff is applied on",
                                                                          "the following selected human populations:")),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('AFKG', 'All MAF KG', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('AMR_AFKG', 'AMR MAF KG', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('ASN_AFKG', 'ASN MAF KG', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('AFR_AFKG', 'AFR MAF KG', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('EUR_AFKG', 'EUR MAF KG', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('AFESP', 'All MAF ESP', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('EA_AFESP', 'EA MAF ESP', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('AA_AFESP', 'AA MAF ESP', TRUE)),
                ## conservation tab
                conditionalPanel(condition="input.tsp == 'conservation'", tags$hr()),
                ## if (!is.na(mtPhastConsDb))
                conditionalPanel(condition="input.tsp == 'conservation'",
                                 checkboxInput('minPhastConsFlag',
                                               'Filter by phastCons score', FALSE)),
                conditionalPanel(condition="input.tsp == 'conservation' && input.minPhastConsFlag == true",
                                 numericInput('minPhastCons', 'Minimum phastCons score:', 0.00)),
                ## if (!is.na(mtGenePhylostrataDb))
                conditionalPanel(condition="input.tsp == 'conservation'",
                                 checkboxInput('minPhylostratumFlag', 'Filter by conserved gene phylostratum', FALSE)),
                conditionalPanel(condition="input.tsp == 'conservation' && input.minPhylostratumFlag == true",
                                 selectInput("minPhylostratum", "Minimum conserved gene phylostratum:",
                                             choices=phylostrata)),
                ## cryptic splice sites tab
                conditionalPanel(condition="input.tsp == 'cryp'", tags$hr()),
                conditionalPanel(condition="input.tsp == 'cryp'", checkboxInput('minCRYP5ssFlag', 'Filter by cryptic 5\'ss', FALSE)),
                conditionalPanel(condition="input.tsp == 'cryp' && input.minCRYP5ssFlag == true",
                                 numericInput('minCRYP5ss', 'Minimum cryptic 5\'ss score:', -99)),
                conditionalPanel(condition="input.tsp == 'cryp'", checkboxInput('minCRYP3ssFlag', 'Filter by cryptic 3\'ss', FALSE)),
                conditionalPanel(condition="input.tsp == 'cryp' && input.minCRYP3ssFlag == true",
                                 numericInput('minCRYP3ss', 'Minimum cryptic 3\'ss score:', -99)),
                tags$hr(),
                downloadButton('downloadData', 'Download Variants'),
                downloadButton('generateReport', 'Generate Report'),
                actionButton('closesavebutton', 'Save & Close')
              ),
              mainPanel(
                tabsetPanel(
                  tabPanel("Genome", tableOutput('tableGenome'), value="genome"),
                  tabPanel("Gene", tableOutput('tableGene'), value="gene"),
                  tabPanel("Transcript", tableOutput('tableTranscript'), value="transcript"),
                  tabPanel("Protein", htmlOutput('tableProtein'), value="protein"),
                  ## if (!is.na(mtMafDb))
                    tabPanel("MAF", tableOutput('tableMAF'), value="maf"),
                  ## if (!is.na(mtPhastConsDb) || !is.na(mtGenePhylostrataDb))
                    tabPanel("Conservation", tableOutput('tableConservation'), value="conservation"),
                  tabPanel("CrypSplice", tableOutput('tableCrypSplice'), value="cryp"),
                  tabPanelAbout(),
                  id="tsp"
                )
              )
            )

  app$server <- function(input, output) {

    observe({
      if (input$closesavebutton == 0)
        return()
      isolate({
        ## individual selected
          if (input$individualFlag)
              selectIndividual(vfResultsObj) <- input$selectIndividual
          else
              selectIndividual(vfResultsObj) <- NA_character_
        ## presence in dbSNP
        if (input$dbSNPpresentFlag)
          dbSNPpresent(vfResultsObj) <- input$dbSNPpresent
        else
          dbSNPpresent(vfResultsObj) <- NA_character_

        ## gene selected
          if (input$geneFlag)
              selectGene(vfResultsObj) <- input$selectGene
          else
              selectGene(vfResultsObj) <- NA_character_
          
        ## presence in OMIM
        if (input$OMIMpresentFlag)
          OMIMpresent(vfResultsObj) <- input$OMIMpresent
        else
          OMIMpresent(vfResultsObj) <- NA_character_

        ## type of variant
        variantType(vfResultsObj) <- input$variantType

        ## type of amino acid change
        aaChangeType(vfResultsObj) <- input$aaChangeType

        ## minimum allele frequency
        if (!is.na(mtMafDb)) {
          mafMask <- MAFpop(vfResultsObj)
          for (i in names(mafMask)) ## unlisting a reactivevalues object does not work :(
            mafMask[i] <- input[[i]]
          MAFpop(vfResultsObj) <- mafMask
          maxMAF(vfResultsObj) <- as.numeric(input$maxMAF)
          naMAF(vfResultsObj) <- input$naMAF
        }

        ## nucleotide conservation
        if (!is.na(mtPhastConsDb)) {
          if (input$minPhastConsFlag)
            minPhastCons(vfResultsObj) <- input$minPhastCons
          else
            minPhastCons(vfResultsObj) <- NA_real_
        }

        ## gene conservation
        if (!is.na(mtGenePhylostrataDb)) {
          if (input$minPhylostratumFlag)
            minPhylostratum(vfResultsObj) <- input$minPhylostratum
          else
            minPhylostratum(vfResultsObj) <- NA_integer_
        }

        ## cryptic splice sites
        if (input$minCRYP5ssFlag)
          minCRYP5ss(vfResultsObj) <- input$minCRYP5ss
        else
          minCRYP5ss(vfResultsObj) <- NA_real_

        if (input$minCRYP3ssFlag)
          minCRYP3ss(vfResultsObj) <- input$minCRYP3ss
        else
          minCRYP3ss(vfResultsObj) <- NA_real_

        stopApp(returnValue=vfResultsObj)
      })
    })

    filteredVariantsReact <- reactive({
        ## individual selected
          if (input$individualFlag)
              selectIndividual(vfResultsObj) <- input$selectIndividual
          else
              selectIndividual(vfResultsObj) <- NA_character_        
      ## presence in dbSNP
      if (input$dbSNPpresentFlag)
        dbSNPpresent(vfResultsObj) <- input$dbSNPpresent
      else
        dbSNPpresent(vfResultsObj) <- NA_character_
          
      ## gene selected
          if (input$geneFlag)
              selectGene(vfResultsObj) <- input$selectGene
          else
              selectGene(vfResultsObj) <- NA_character_
          
      ## presence in OMIM
      if (input$OMIMpresentFlag)
        OMIMpresent(vfResultsObj) <- input$OMIMpresent
      else
        OMIMpresent(vfResultsObj) <- NA_character_

      ## type of variant
      variantType(vfResultsObj) <- input$variantType

      ## type of amino acid change
      aaChangeType(vfResultsObj) <- input$aaChangeType

      ## minimum allele frequency
      if (!is.na(mtMafDb)) {
        mafMask <- MAFpop(vfResultsObj)
        for (i in names(mafMask)) ## unlisting a reactivevalues object does not work :(
          mafMask[i] <- input[[i]]
        MAFpop(vfResultsObj) <- mafMask
        maxMAF(vfResultsObj) <- as.numeric(input$maxMAF)
        naMAF(vfResultsObj) <- input$naMAF
      }

      ## nucleotide conservation
      if (!is.na(mtPhastConsDb)) {
        if (input$minPhastConsFlag)
          minPhastCons(vfResultsObj) <- input$minPhastCons
        else
          minPhastCons(vfResultsObj) <- NA_real_
      }

      ## gene conservation
      if (!is.na(mtGenePhylostrataDb)) {
        if (input$minPhylostratumFlag)
          minPhylostratum(vfResultsObj) <- input$minPhylostratum
        else
          minPhylostratum(vfResultsObj) <- NA_integer_
      }

      ## cryptic splice sites
      if (input$minCRYP5ssFlag)
        minCRYP5ss(vfResultsObj) <- as.numeric(input$minCRYP5ss)
      else
        minCRYP5ss(vfResultsObj) <- NA_real_

      if (input$minCRYP3ssFlag)
        minCRYP3ss(vfResultsObj) <- as.numeric(input$minCRYP3ss)
      else
        minCRYP3ss(vfResultsObj) <- NA_real_

      vdf <- filteredVariants(vfResultsObj)
      vdf <- DataFrame(VarID=names(vdf), CHR=seqnames(vdf),
                       POS=start(vdf), POSITION=start(vdf),
                       mcols(vdf))
      ## the ALT column is a DNAStringSetList object and requires this
      ## step to coerce to character
      vdf$ALT <- sapply(vdf$ALT, paste, collapse=",")
      vdf <- as.data.frame(vdf)

      ## both Homozygous and Heterozygous columns are CharacterList objects and
      ## requires this step to coerce to character
      vdf$HomozygousAlt <- sapply(vdf$HomozygousAlt, paste, collapse=", ")
      vdf$Heterozygous <- sapply(vdf$Heterozygous, paste, collapse=", ")
      
      if (nrow(vdf) > 0) {
        varlocs <- paste0("<a href=http://genome.ucsc.edu/cgi-bin/hgTracks?org=", species,
                          "&db=", genomeVersion, 
                          "&position=", vdf$CHR, ":", vdf$POS, " target=\"ucsc\">",
                          vdf$CHR, ":", vdf$POS,
                          "</a>")
        tempOMIM <- rep(NA_character_, length(vfResultsObj))
        tempOMIM[!is.na(vdf$OMIM)] <- sapply(strsplit(vdf$OMIM[!is.na(vdf$OMIM)], ","), function(x) {
          z <- paste0("<a href=http://www.omim.org/entry/", x, " target=\"omim\">", x, "</a>")
          a <- paste(z, collapse=", ")
          a
        })
      
        tempdbSNP <- rep(NA_character_, length(vfResultsObj))
        nors <- sub("rs", "", vdf$dbSNP[!is.na(vdf$dbSNP)])
        tempdbSNP[!is.na(vdf$dbSNP)] <- paste0("<a href=http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=",
                           nors, " target=\"dbsnp\">", vdf$dbSNP[!is.na(vdf$dbSNP)], "</a>")
      
        vdf[["POSITION"]] <- varlocs
        vdf[["OMIM"]] <- tempOMIM
        vdf[["dbSNP"]] <- tempdbSNP
      }

      vdf <- vdf[, -match(c("CHR", "POS"), colnames(vdf))]
      
      rownames(vdf) <- NULL

      vdf 
    })

    output$tableGenome <- renderTable({
      filteredVariantsReact()[, c("VarID", "POSITION", "dbSNP", "TYPE", "Heterozygous", "HomozygousAlt")]
    }, NA.string="NA",  sanitize.text.function=function(x){x})

    output$tableGene <- renderTable({
      filteredVariantsReact()[, c("VarID", "POSITION", "GENE", "LOCATION", "OMIM")]
    }, NA.string="NA",  sanitize.text.function=function(x){x})

    output$tableTranscript <- renderTable({
      filteredVariantsReact()[, c("VarID", "POSITION", "GENE", "TXID", "LOCATION", "LOCSTART", "cDNALOC", "CDS")]
    }, NA.string="NA",  sanitize.text.function=function(x){x})

    output$tableProtein <- renderTable({
      filteredVariantsReact()[, c("VarID", "GENE", "CONSEQUENCE", "AAchange", "AAchangeType", "PolyPhen2", "PROVEAN")]
    }, NA.string="NA",  sanitize.text.function=function(x){x})

    output$tableMAF <- renderTable({
      selcnames <- c("VarID", "dbSNP", "Max", "All KG", "AMR KG", "ASN KG",
                     "AFR KG", "EUR KG", "All ESP", "EA ESP", "AA ESP")
      fv <- filteredVariantsReact()
      cnamesmt <- match(c("VarID", "dbSNP", "maxMAF", "AFKG", "AMR_AFKG", "ASN_AFKG",
                          "AFR_AFKG", "EUR_AFKG", "AFESP", "EA_AFESP", "AA_AFESP"),
                        colnames(fv))
      fv <- fv[, na.omit(cnamesmt)]
      colnames(fv) <- selcnames[!is.na(cnamesmt)]
      fv
    }, NA.string="NA",  sanitize.text.function=function(x){x})

    output$tableConservation <- renderTable({
      fv <- filteredVariantsReact()[, c("VarID", "POSITION", "GENE", "phastCons", "GenePhylostratum", "GenePhylostratumTaxID")]
      fv$GenePhylostratum[!is.na(fv$GenePhylostratum)] <- sprintf("<a href=\"http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s&lvl=5\" target='ncbitaxonomybrowser'>%s</a>", fv$GenePhylostratumTaxID[!is.na(fv$GenePhylostratum)], fv$GenePhylostratum[!is.na(fv$GenePhylostratum)])
      fv <- fv[, -match("GenePhylostratumTaxID", colnames(fv))]
      fv
    }, NA.string="NA", sanitize.text.function=function(x){x})

    output$tableCrypSplice <- renderTable({
      colscrypss <- c("CRYP5ssREF", "CRYP5ssALT", "CRYP5ssPOS", "CRYP3ssREF", "CRYP3ssALT", "CRYP3ssPOS")
      selcnames <- c("VarID", "POSITION", "5'ss Ref", "5'ss Alt", "5'ss Pos", "3'ss Ref", "3'ss Alt", "3'ss Pos")

      fv <- filteredVariantsReact()[, c("VarID", "POSITION", colscrypss), drop=FALSE]
      colnames(fv) <- selcnames
      fv
    }, NA.string="NA",  sanitize.text.function=function(x){x})

    ## this provides a tab-separated value file with the selected vfResultsObj
    output$downloadData <- downloadHandler(
      filename = function() { "vfResultsObj.tsv" },
      content = function(file) {
        ## individual selected
          if (input$individualFlag)
              selectIndividual(vfResultsObj) <- input$selectIndividual
          else
              selectIndividual(vfResultsObj) <- NA_character_        
        
        ## presence in dbSNP
        if (input$dbSNPpresentFlag)
          dbSNPpresent(vfResultsObj) <- input$dbSNPpresent
        else
          dbSNPpresent(vfResultsObj) <- NA_character_

        ## gene selected
          if (input$geneFlag)
              selectGene(vfResultsObj) <- input$selectGene
          else
              selectGene(vfResultsObj) <- NA_character_
          
        ## presence in OMIM
        if (input$OMIMpresentFlag)
          OMIMpresent(vfResultsObj) <- input$OMIMpresent
        else
          OMIMpresent(vfResultsObj) <- NA_character_

        ## type of variant
        variantType(vfResultsObj) <- input$variantType

        ## type of amino acid change
        aaChangeType(vfResultsObj) <- input$aaChangeType

        ## minimum allele frequency
        if (!is.na(mtMafDb)) {
          mafMask <- MAFpop(vfResultsObj)
          for (i in names(mafMask)) ## unlisting a reactivevalues object does not work :(
            mafMask[i] <- input[[i]]
          MAFpop(vfResultsObj) <- mafMask
          maxMAF(vfResultsObj) <- as.numeric(input$maxMAF)
          naMAF(vfResultsObj) <- input$naMAF
        }

        ## nucleotide conservation
        if (!is.na(mtPhastConsDb)) {
          if (input$minPhastConsFlag)
            minPhastCons(vfResultsObj) <- input$minPhastCons
          else
            minPhastCons(vfResultsObj) <- NA_real_
        }

        ## gene conservation
        if (!is.na(mtGenePhylostrataDb)) {
          if (input$minPhylostratumFlag)
            minPhylostratum(vfResultsObj) <- input$minPhylostratum
          else
            minPhylostratum(vfResultsObj) <- NA_integer_
        }

        ## cryptic splice sites
        if (input$minCRYP5ssFlag)
          minCRYP5ss(vfResultsObj) <- as.numeric(input$minCRYP5ss)
        else
          minCRYP5ss(vfResultsObj) <- NA_real_

        if (input$minCRYP3ssFlag)
          minCRYP3ss(vfResultsObj) <- as.numeric(input$minCRYP3ss)
        else
          minCRYP3ss(vfResultsObj) <- NA_real_

        vdf <- filteredVariants(vfResultsObj)
        vdf <- DataFrame(VarID=names(vdf), CHR=seqnames(vdf),
                         POS=start(vdf), POSITION=start(vdf),
                         mcols(vdf))
        ## the ALT column is a DNAStringSetList object and requires this
        ## step to coerce to character
        vdf$ALT <- sapply(vdf$ALT, paste, collapse=",")
        vdf <- as.data.frame(vdf)
      
        write.table(vdf, file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
      },
      contentType = "text/plain"
    )

    ## this should generate a report to facilitate reproducing the analysis with R code
    output$generateReport <- downloadHandler(
      filename = function() { "report.txt" },
      content = function(file) {
        ## individual selected
          if (input$individualFlag)
              selectIndividual(vfResultsObj) <- input$selectIndividual
          else
              selectIndividual(vfResultsObj) <- NA_character_        

        ## presence in dbSNP
        if (input$dbSNPpresentFlag)
          dbSNPpresent(vfResultsObj) <- input$dbSNPpresent
        else
          dbSNPpresent(vfResultsObj) <- NA_character_

        ## gene selected
          if (input$geneFlag)
              selectGene(vfResultsObj) <- input$selectGene
          else
              selectGene(vfResultsObj) <- NA_character_
          
        ## presence in OMIM
        if (input$OMIMpresentFlag)
          OMIMpresent(vfResultsObj) <- input$OMIMpresent
        else
          OMIMpresent(vfResultsObj) <- NA_character_

        ## type of variant
        variantType(vfResultsObj) <- input$variantType

        ## type of amino acid change
        aaChangeType(vfResultsObj) <- input$aaChangeType

        ## minimum allele frequency
        if (!is.na(mtMafDb)) {
          mafMask <- MAFpop(vfResultsObj)
          for (i in names(mafMask)) ## unlisting a reactivevalues object does not work :(
            mafMask[i] <- input[[i]]
          MAFpop(vfResultsObj) <- mafMask
          maxMAF(vfResultsObj) <- as.numeric(input$maxMAF)
          naMAF(vfResultsObj) <- input$naMAF
        }

        ## nucleotide conservation
        if (!is.na(mtPhastConsDb)) {
          if (input$minPhastConsFlag)
            minPhastCons(vfResultsObj) <- input$minPhastCons
          else
            minPhastCons(vfResultsObj) <- NA_real_
        }

        ## gene conservation
        if (!is.na(mtGenePhylostrataDb)) {
          if (input$minPhylostratumFlag)
            minPhylostratum(vfResultsObj) <- input$minPhylostratum
          else
            minPhylostratum(vfResultsObj) <- NA_integer_
        }

        ## cryptic splice sites
        if (input$minCRYP5ssFlag)
          minCRYP5ss(vfResultsObj) <- as.numeric(input$minCRYP5ss)
        else
          minCRYP5ss(vfResultsObj) <- NA_real_

        if (input$minCRYP3ssFlag)
          minCRYP3ss(vfResultsObj) <- as.numeric(input$minCRYP3ss)
        else
          minCRYP3ss(vfResultsObj) <- NA_real_


        sink(file)
        cat("R script and output to get the resulting filtered variants with VariantFiltering\n")
        cat("================================================================================\n\n")

        cat("> library(VariantFiltering)\n\n")
        cat(sprintf("> %s <- %s\n", gettext(vfResultsObj@callObj)[2], param(vfResultsObj)$callStr))          
        cat(sprintf("> res <- %s\n", deparse(vfResultsObj@callObj)))
        if (!is.na(dbSNPpresent(vfResultsObj)))
          cat(sprintf("> dbSNPpresent(res) <- \"%s\"\n", dbSNPpresent(vfResultsObj)))
        if (!is.na(OMIMpresent(vfResultsObj)))
          cat(sprintf("> OMIMpresent(res) <- \"%s\"\n", OMIMpresent(vfResultsObj)))
        if (variantType(vfResultsObj) != "Any")
          cat(sprintf("> variantType(res) <- \"%s\"\n", variantType(vfResultsObj)))
        if (aaChangeType(vfResultsObj) != "Any")
          cat(sprintf("> aaChangeType(res) <- \"%s\"\n", aaChangeType(vfResultsObj)))
        if (!is.na(mtMafDb)) {
          if (!all(MAFpop(vfResultsObj))) ## default values need not to be set
            cat(sprintf("> MAFpop(res) <- c(%s)\n", paste(paste(names(MAFpop(vfResultsObj)), as.character(MAFpop(vfResultsObj)), sep="="), collapse=", ")))
          if (!naMAF(vfResultsObj)) ## default values need not to be set
            cat(sprintf("> naMAF(res) <- %s\n", ifelse(input$naMAF, "TRUE", "FALSE")))
          cat(sprintf("> maxMAF(res) <- %.3f\n", input$maxMAF))
        }
        if (!is.na(mtPhastConsDb)) {
          if (!is.na(minPhastCons(vfResultsObj)))
            cat(sprintf("> minPhastCons(res) <- %.1f\n", minPhastCons(vfResultsObj)))
        }
        if (!is.na(mtGenePhylostrataDb)) {
          if (!is.na(minPhylostratum(vfResultsObj)))
            cat(sprintf("> minPhylostratum(res) <- %d\n", minPhylostratum(vfResultsObj)))
        }
        if (!is.na(minCRYP5ss(vfResultsObj)))
          cat(sprintf("> minCRYP5ss(res) <- %.2f\n", minCRYP5ss(vfResultsObj)))
        if (!is.na(minCRYP3ss(vfResultsObj)))
          cat(sprintf("> minCRYP3ss(res) <- %.2f\n", minCRYP3ss(vfResultsObj)))
        cat("> res\n")
        print(vfResultsObj)
        cat("> reportVariants(res, type=\"tsv\", file=\"variants.tsv\")\n")
        cat("> sessionInfo()\n")
        print(sessionInfo())
        sink()
      },
      contentType = "text/plain"
    )

  }

  runApp(app)
})
