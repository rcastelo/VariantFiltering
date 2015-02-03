## group by variant position when providing the summary by variant (??????)

setMethod("show", signature(object="VariantFilteringResults"),
          function(object) {
            cat("\nVariantFiltering results object\n\n")
            if (inheritanceModel(object) != "any")
              cat(sprintf("  Variants segregate according to a/an %s inheritance model\n", inheritanceModel(object)))
            else
              cat("  Variants are not filtered by inheritance model\n")
            cat("  Functional annotation filters\n")
            if (is.na(dbSNPpresent(object)))
                cat(sprintf("    No filtering on presence in %s %s\n", provider(param(object)$snpdb),
                            releaseName(param(object)$snpdb)))
            else
                cat(sprintf("    Present in %s %s: %s\n", provider(param(object)$snpdb),
                            releaseName(param(object)$snpdb), dbSNPpresent(object)))
            cat(sprintf("    Variant type: %s\n", variantType(object)))
            cat(sprintf("    Amino acid change type: %s\n", aaChangeType(object)))
            if ("MafDb" %in% sapply(param(object)$otherAnnotations, class)) {
              cat(sprintf("    Populations used for MAF filtering: %s\n", paste(names(MAFpop(object))[MAFpop(object)], collapse=", ")))
              cat(sprintf("    Include MAF NA values: %s\n", ifelse(naMAF(object), "yes", "no")))
              cat(sprintf("    Maximum MAF: %.2f\n", maxMAF(object)))
            }
            if ("PhastConsDb" %in% sapply(param(object)$otherAnnotations, class)) {
              if (is.na(minPhastCons(object)))
                cat("    No filtering on nucleotide conservation\n")
              else
                cat(sprintf("    Minimum score for phastCons nucleotide conservation: %.2f\n", minPhastCons(object)))
            }
            if ("GenePhylostrataDb" %in% sapply(param(object)$otherAnnotations, class)) {
              whGPSdb <- match("GenePhylostrataDb", sapply(param(object)$otherAnnotations, class))
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
            if (!is.na(match("intergenic", names(vxloc))))
              cat(sprintf("    %s (%.1f%%) located in intergenic regions\n", sprintf(paddgts, vxloc["intergenic"]), 100*vxloc["intergenic"]/total))
          })

setMethod("param", signature(x="VariantFilteringResults"),
          function(x) {
            x@inputParameters
          })

setMethod("inheritanceModel", signature(x="VariantFilteringResults"),
          function(x) {
            x@inheritanceModel
          })

setMethod("dbSNPpresent", signature(x="VariantFilteringResults"),
          function(x) {
            x@dbSNPflag
          })

setReplaceMethod("dbSNPpresent", signature(x="VariantFilteringResults", value="ANY"),
                 function(x, value) {
                   if (!is.na(value) && !is.character(value) && !is.logical(value))
                     stop("Presence in dbSNP should be indicated either by NA, a character string \"yes\" or \"no\", or a logical value TRUE or FALSE.")

                   if (!is.na(value)) {
                     if (is.character(value)) {
                       if (tolower(value) == "yes")
                         value <- "Yes"
                       else if (tolower(value) == "no")
                         value <- "No"
                       else
                         stop("Presence in dbSNP should be indicated either by NA, a character string \"yes\" or \"no\", or a logical value TRUE or FALSE.")
                     }

                     if (is.logical(value)) {
                       if (value)
                         value <- "Yes"
                       else
                         value <- "No"
                     }
                   }

                   x@dbSNPflag <- as.character(value)

                   x
                 })

setMethod("OMIMpresent", signature(x="VariantFilteringResults"),
          function(x) {
            x@OMIMflag
          })

setReplaceMethod("OMIMpresent", signature(x="VariantFilteringResults", value="ANY"),
                 function(x, value) {
                   if (!is.na(value) && !is.character(value) && !is.logical(value))
                     stop("Presence in OMIM should be indicated either by NA, a character string \"yes\" or \"no\", or a logical value TRUE or FALSE.")

                   if (!is.na(value)) {
                     if (is.character(value)) {
                       if (tolower(value) == "yes")
                         value <- "Yes"
                       else if (tolower(value) == "no")
                         value <- "No"
                       else
                         stop("Presence in OMIM should be indicated either by NA, a character string \"yes\" or \"no\", or a logical value TRUE or FALSE.")
                     }

                     if (is.logical(value)) {
                       if (value)
                         value <- "Yes"
                       else
                         value <- "No"
                     }
                   }

                   value <- as.character(value)
                   x@OMIMflag <- value

                   x
                 })

setMethod("variantType", signature(x="VariantFilteringResults"),
          function(x) {
            x@variantType
          })

setReplaceMethod("variantType", signature(x="VariantFilteringResults", value="character"),
                 function(x, value=c("Any", "SNV", "InDel", "MNV")) {
                   value <- match.arg(value)
                   x@variantType <- value
                   x
                 })

setMethod("aaChangeType", signature(x="VariantFilteringResults"),
          function(x) {
            x@aaChangeType
          })

setReplaceMethod("aaChangeType", signature(x="VariantFilteringResults", value="character"),
                 function(x, value=c("Any", "Radical", "Conservative")) {
                   value <- match.arg(value)
                   x@aaChangeType <- value
                   x
                 })
setMethod("naMAF", signature(x="VariantFilteringResults"),
          function(x) {
            if (is.na(match("MafDb", sapply(param(x)$otherAnnotations, class))))
              stop("A MafDb object was not used to annotate variants.")

            x@naMAF
          })

setReplaceMethod("naMAF", signature(x="VariantFilteringResults", value="logical"),
          function(x, value) {
            if (is.na(match("MafDb", sapply(param(x)$otherAnnotations, class))))
              stop("A MafDb object was not used to annotate variants.")

            x@naMAF <- value
            x
          })

setMethod("MAFpop", signature(x="VariantFilteringResults"),
          function(x) {
            if (is.na(match("MafDb", sapply(param(x)$otherAnnotations, class))))
              stop("A MafDb object was not used to annotate variants.")

            x@MAFpopMask
          })

setReplaceMethod("MAFpop", signature(x="VariantFilteringResults", value="logical"),
                 function(x, popkey=NA, value) {
                   if (is.na(match("MafDb", sapply(param(x)$otherAnnotations, class))))
                     stop("A MafDb object was not used to annotate variants.")

                   if (is.na(popkey)) {
                     if (is.null(names(value)))
                       stop("The given value must be a logical vector whose names match the available population keywords.")
                     if (any(is.na(match(names(value), names(x@MAFpopMask)))))
                       stop(sprintf("Element names %s do not match the available population keywords",
                                    names(value)[is.na(match(names(value), names(x@MAFpopMask)))]))
                     popkey <- names(value)
                   } else {
                     if (any(is.na(match(popkey, names(x@MAFpopMask)))))
                       stop(sprintf("%s does not match the available population keywords", popkey))
                   }

                   x@MAFpopMask[popkey] <- value
                   x
                 })

setMethod("maxMAF", signature(x="VariantFilteringResults"),
          function(x) {
            if (is.na(match("MafDb", sapply(param(x)$otherAnnotations, class))))
              stop("A MafDb object was not used to annotate variants.")

            x@maxMAF
          })

setReplaceMethod("maxMAF", signature(x="VariantFilteringResults", value="numeric"),
                 function(x, value) {
                   if (is.na(match("MafDb", sapply(param(x)$otherAnnotations, class))))
                     stop("A MafDb object was not used to annotate variants.")

                   x@maxMAF <- value
                   x
                 })

setMethod("minPhastCons", signature(x="VariantFilteringResults"),
          function(x) {
            if (is.na(match("PhastConsDb", sapply(param(x)$otherAnnotations, class))))
              stop("A phastConsDb object was not used to annotate variants.")

            x@minPhastCons
          })

setReplaceMethod("minPhastCons", signature(x="VariantFilteringResults", value="ANY"),
                 function(x, value) {
                   if (is.na(match("PhastConsDb", sapply(param(x)$otherAnnotations, class))))
                     stop("A phastConsDb object was not used to annotate variants.")

                   if (!is.na(value) && !is.numeric(value) && !is.integer(value))
                     stop("Only a numeric or NA value is allowed as minimum cutoff for phastCons scores.")

                   if (!is.na(value)) {
                     if (value < 0 || value > 1)
                       stop("The minimum cutoff for phastCons scores should be a number between 0 and 1.")
                   }

                   value <- as.numeric(value)
                   x@minPhastCons <- value

                   x
                 })

setMethod("minPhylostratum", signature(x="VariantFilteringResults"),
          function(x) {
            if (is.na(match("GenePhylostrataDb", sapply(param(x)$otherAnnotations, class))))
              stop("A GenePhylostrataDb object was not used to annotate variants.")

            x@minPhylostratumIndex
          })

setReplaceMethod("minPhylostratum", signature(x="VariantFilteringResults", value="ANY"),
                 function(x, value) {
                   whGPSdb <- match("GenePhylostrataDb", sapply(param(x)$otherAnnotations, class))
                   if (is.na(whGPSdb))
                     stop("A GenePhylostrataDb object was not used to annotate variants.")

                   if (!is.na(value) && !is.numeric(value) && !is.integer(value) && !is.character(value))
                     stop("Only an integer, character or NA value is allowed as minimum phylostratum index.")

                   if (!is.na(value)) {
                     if (is.character(value)) {
                       value <- match(tolower(value),
                                      tolower(genePhylostrata(param(x)$otherAnnotations[[whGPSdb]])$Description))
                       if (is.na(value))
                         stop(sprintf("%s is not valid. The minimum phylostratum character value should be one of: %s",
                                      value, paste(genePhylostrata(param(x)$otherAnnotations[[whGPSdb]])$Description,
                                                   collapse=",")))
                     }

                     value <- as.integer(value)
                     if (value < 1 || value > nrow(genePhylostrata(param(x)$otherAnnotations[[whGPSdb]])))
                       stop(sprintf("%d is not valid. The minimum phylostratum integer value should be between 1 and %d",
                                    nrow(genePhylostrata(param(x)$otherAnnotations[[whGPSdb]]))))
                   }

                   value <- as.integer(value)
                   x@minPhylostratumIndex <- value

                   x
                 })

setMethod("minCRYP5ss", signature(x="VariantFilteringResults"),
          function(x) {
            x@minCRYP5ss
          })

setReplaceMethod("minCRYP5ss", signature(x="VariantFilteringResults", value="ANY"),
                 function(x, value) {
                   if (!is.na(value) && !is.numeric(value) && !is.integer(value))
                     stop("Only a numeric or NA value is allowed as minimum cutoff for cryptic 5'ss.")
                   x@minCRYP5ss <- as.numeric(value)
                   x
                 })

setMethod("minCRYP3ss", signature(x="VariantFilteringResults"),
          function(x) {
            x@minCRYP3ss
          })

setReplaceMethod("minCRYP3ss", signature(x="VariantFilteringResults", value="ANY"),
                 function(x, value) {
                   if (!is.na(value) && !is.numeric(value) && !is.integer(value))
                     stop("Only a numeric or NA value is allowed as minimum cutoff for cryptic 3'ss.")
                   x@minCRYP3ss <- as.numeric(value)
                   x
                 })


## get all variants without applying any filter
setMethod("allVariants", signature(x="VariantFilteringResults"), 
          function(x) {
            x@variants
          })

## get variants after applying all filters
setMethod("filteredVariants", signature(x="VariantFilteringResults"), 
          function(x, unusedColumns.rm=FALSE) {
            vars <- allVariants(x)
            rowsMask <- rep(TRUE, length(vars))

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

            vars[rowsMask, colsMask, drop=FALSE]
          })

## shiny app to filter and visualize variants
setMethod("reportVariants", signature(vfResultsObj="VariantFilteringResults"),
          function(vfResultsObj, type=c("shiny", "csv", "tsv"), file=NULL) {
              
  type <- match.arg(type)

  if (class(vfResultsObj) != "VariantFilteringResults")
    stop("Input argument 'vfResultsObj' should be an object of class 'VariantFilteringResults'.")

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
  
  ## parameters for the UCSC genome browser. the organism name is derived from
  ## the 'organism()' getter for BSgenome objects and 
  org <- organism(param(vfResultsObj)$bsgenome)
  org <- paste0(substr(org, 1, 1), strsplit(org, " ")[[1]][2])
  genomeBuild <- providerVersion(param(vfResultsObj)$bsgenome)

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
                                 checkboxInput('dbSNPpresentFlag',
                                               'Filter by presence in dbSNP', FALSE)),
                conditionalPanel(condition="input.tsp == 'genome' && input.dbSNPpresentFlag == true",
                                 selectInput("dbSNPpresent", "Present in dbSNP:",
                                             choices=c("Yes", "No"))),
                conditionalPanel(condition="input.tsp == 'genome'", selectInput("variantType", "Variant Type:",
                                                                    choices=c("Any", "SNV", "InDel", "MNV"))),
                ## gene tab
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
        ## presence in dbSNP
        if (input$dbSNPpresentFlag)
          dbSNPpresent(vfResultsObj) <- input$dbSNPpresent
        else
          dbSNPpresent(vfResultsObj) <- NA_character_

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
      ## presence in dbSNP
      if (input$dbSNPpresentFlag)
        dbSNPpresent(vfResultsObj) <- input$dbSNPpresent
      else
        dbSNPpresent(vfResultsObj) <- NA_character_

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
      vdf <- DataFrame(VarID=vdf$VARID, CHR=seqnames(vdf),
                       POS=start(vdf), POSITION=start(vdf),
                       mcols(vdf))
      ## the ALT column is a DNAStringSetList object and requires this
      ## step to coerce to character
      ## vdf$ALT <- sapply(vdf$ALT, paste, collapse=",")
      vdf$REF <- ref(vdf)
      vdf$ALT <- alt(vdf)
      vdf <- as.data.frame(vdf)
      
      if (nrow(vdf) > 0) {
        varlocs <- paste0("<a href=http://genome.ucsc.edu/cgi-bin/hgTracks?org=", species,
                          "&db=", genomeBuild, 
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
      filteredVariantsReact()[, c("VarID", "POSITION", "dbSNP", "TYPE")]
    }, NA.string="NA",  sanitize.text.function=function(x){x})

    output$tableGene <- renderTable({
      filteredVariantsReact()[, c("VarID", "POSITION", "GENE", "LOCATION", "OMIM")]
    }, NA.string="NA",  sanitize.text.function=function(x){x})

    output$tableTranscript <- renderTable({
      ## filteredVariantsReact()[, c("VarID", "POSITION", "GENE", "TXID", "LOCATION", "LOCSTART", "cDNALOC", "CDS")]
      filteredVariantsReact()[, c("VarID", "POSITION", "GENE", "TXID", "LOCATION", "LOCSTART", "cDNALOC", "CDS")]
    }, NA.string="NA",  sanitize.text.function=function(x){x})

    output$tableProtein <- renderTable({
      ## filteredVariantsReact()[, c("VarID", "GENE", "CONSEQUENCE", "AAchange", "AAchangeType", "PolyPhen2", "PROVEAN")]
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
        ## presence in dbSNP
        if (input$dbSNPpresentFlag)
          dbSNPpresent(vfResultsObj) <- input$dbSNPpresent
        else
          dbSNPpresent(vfResultsObj) <- NA_character_

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
          maxMAF(vfResultsObj) <- input$maxMAF
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
        ## presence in dbSNP
        if (input$dbSNPpresentFlag)
          dbSNPpresent(vfResultsObj) <- input$dbSNPpresent
        else
          dbSNPpresent(vfResultsObj) <- NA_character_

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
          maxMAF(vfResultsObj) <- input$maxMAF
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
