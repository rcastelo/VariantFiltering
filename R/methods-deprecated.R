##
## deprecated methods, to become defunct in BioC 3.8
##

setMethod("dbSNPpresent", signature(x="VariantFilteringResults"),
          function(x) {
            warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
            active(filters(x))["dbSNP"]
          })

setReplaceMethod("dbSNPpresent", signature(x="VariantFilteringResults", value="ANY"),
                 function(x, value) {
                   warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
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

                   ## x@dbSNPflag <- as.character(value)
                   if (value == "Yes")
                     active(filters(x))["dbSNP"] <- TRUE
                   else
                     active(filters(x))["dbSNP"] <- FALSE

                   x
                 })

setMethod("OMIMpresent", signature(x="VariantFilteringResults"),
          function(x) {
            warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
            active(filters(x))["OMIM"]
          })

setReplaceMethod("OMIMpresent", signature(x="VariantFilteringResults", value="ANY"),
                 function(x, value) {
                   warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")

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

                   ## value <- as.character(value)
                   ## x@OMIMflag <- value
                   if (value == "Yes")
                     active(filters(x))["OMIM"] <- TRUE
                   else
                     active(filters(x))["OMIM"] <- FALSE

                   x
                 })

setMethod("variantType", signature(x="VariantFilteringResults"),
          function(x) {
            warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
            cutoffs(x)[["variantType"]]
          })

setReplaceMethod("variantType", signature(x="VariantFilteringResults", value="logical"),
                 function(x, typkey=NA, value) {
                   warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")

                   if (any(is.na(value)))
                     stop("The given value(s) must be either TRUE or FALSE.")

                   if (is.na(typkey)) {
                     typkey <- names(cutoffs(x)[["variantType"]])
                     if (is.null(names(value))) {
                       if (length(value) > 1)
                         stop("When given multiple values they must be a logical vector whose names match the available location keywords.")
                       value <- do.call("names<-", list(rep(value, length(typkey)), typkey))
                     } else if (any(is.na(match(names(value), names(cutoffs(x)[["variantType"]])))))
                       stop(sprintf("Element names %s do not match the available variant type keywords",
                                    names(value)[is.na(match(names(value), names(cutoffs(x)[["variantType"]])))]))
                     else
                       typkey <- names(value)
                   } else {
                     if (any(is.na(match(typkey, names(cutoffs(x)[["variantType"]])))))
                       stop(sprintf("%s does not match the available variant type keywords", typkey))
                   }

                   ## x@variantTypeMask[typkey] <- value
                   change(cutoffs(x), "variantType") <- value
                   x
                 })

setMethod("aaChangeType", signature(x="VariantFilteringResults"),
          function(x) {
            warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
            cutoffs(x)[["aaChangeType"]]
          })

setReplaceMethod("aaChangeType", signature(x="VariantFilteringResults", value="character"),
                 function(x, value=c("Any", "Radical", "Conservative")) {
                   warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
                   value <- match.arg(value)
                   ## x@aaChangeType <- value
                   change(cutoffs(x), "aaChangeType") <- value
                   x
                 })

setMethod("variantLocation", signature(x="VariantFilteringResults"),
          function(x) {
            warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
            cutoffs(x)[["location"]]
          })

setReplaceMethod("variantLocation", signature(x="VariantFilteringResults", value="logical"),
                 function(x, lockey=NA, value) {
                   warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
                   if (any(is.na(value)))
                     stop("The given value(s) must be either TRUE or FALSE.")

                   if (is.na(lockey)) {
                     lockey <- names(cutoffs(x)[["location"]])
                     if (is.null(names(value))) {
                       if (length(value) > 1)
                         stop("When given multiple values they must be a logical vector whose names match the available location keywords.")
                       value <- do.call("names<-", list(rep(value, length(lockey)), lockey))
                     } else if (any(is.na(match(names(value), names(cutoffs(x)[["location"]])))))
                       stop(sprintf("Element names %s do not match the available location keywords",
                                    names(value)[is.na(match(names(value), names(cutoffs(x)[["location"]])))]))
                     else
                       lockey <- names(value)
                   } else {
                     if (any(is.na(match(lockey, names(cutoffs(x)[["location"]])))))
                       stop(sprintf("%s does not match the available location keywords", lockey))
                   }

                   ## x@locationMask[lockey] <- value
                   change(cutoffs(x), "location") <- value
                   x
                 })

setMethod("variantConsequence", signature(x="VariantFilteringResults"),
          function(x) {
            warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
            cutoffs(x)[["consequence"]]
          })

setReplaceMethod("variantConsequence", signature(x="VariantFilteringResults", value="logical"),
                 function(x, conkey=NA, value) {
                   warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")

                   if (any(is.na(value)))
                     stop("The given value(s) must be either TRUE or FALSE.")

                   if (is.na(conkey)) {
                     conkey <- names(cutoffs(x)[["consequence"]])
                     if (is.null(names(value))) {
                       if (length(value) > 1)
                         stop("When given multiple values they must be a logical vector whose names match the available consequence keywords.")
                       value <- do.call("names<-", list(rep(value, length(conkey)), conkey))
                     } else if (any(is.na(match(names(value), names(cutoffs(x)[["consequence"]])))))
                       stop(sprintf("Element names %s do not match the available consequence keywords",
                                    names(value)[is.na(match(names(value), names(cutoffs(x)[["consequence"]])))]))
                     else
                       conkey <- names(value)
                   } else {
                     if (any(is.na(match(conkey, names(cutoffs(x)[["consequence"]])))))
                       stop(sprintf("%s does not match the available consequence keywords", conkey))
                   }

                   ## x@consequenceMask[conkey] <- value
                   change(cutoffs(x), "consequence") <- value
                   x
                 })

setMethod("naMAF", signature(x="VariantFilteringResults"),
          function(x) {
            warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
            if (!"MafDb" %in% param(x)$otherAnnotationsClass)
              stop("A MafDb object was not used to annotate variants.")

            active(filters(x))["maxMAF"]
          })

setReplaceMethod("naMAF", signature(x="VariantFilteringResults", value="logical"),
                 function(x, value) {
                   warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
                   if (!"MafDb" %in% param(x)$otherAnnotationsClass)
                     stop("A MafDb object was not used to annotate variants.")

                   ## x@naMAF <- value
                   active(filters(x))["maxMAF"] <- value
                   x
                 })

setMethod("MAFpop", signature(x="VariantFilteringResults"),
          function(x) {
            warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
            if (!"MafDb" %in% param(x)$otherAnnotationsClass)
              stop("A MafDb object was not used to annotate variants.")

            cutoffs(x)[["maxMAF"]]$popmask
          })

setReplaceMethod("MAFpop", signature(x="VariantFilteringResults", value="logical"),
                 function(x, popkey=NA, value) {
                   warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")

                   if (!"MafDb" %in% param(x)$otherAnnotationsClass)
                     stop("A MafDb object was not used to annotate variants.")

                   if (any(is.na(value)))
                     stop("The given value(s) must be either TRUE or FALSE.")

                   if (is.na(popkey)) {
                     popkey <- names(cutoffs(x)[["maxMAF"]]$popmask)
                     if (is.null(names(value))) {
                       if (length(value) > 1)
                         stop("When given multiple values they must be a logical vector whose names match the available population keywords.")
                       value <- do.call("names<-", list(rep(value, length(popkey)), popkey))
                     } else if (any(is.na(match(names(value), names(cutoffs(x)[["maxMAF"]]$popmask)))))
                       stop(sprintf("Element names %s do not match the available population keywords",
                                    names(value)[is.na(match(names(value), names(cutoffs(x)[["maxMAF"]]$popmask)))]))
                     else
                       popkey <- names(value)
                   } else {
                     if (any(is.na(match(popkey, names(cutoffs(x)[["maxMAF"]]$popmask)))))
                       stop(sprintf("%s does not match the available population keywords", popkey))
                   }

                   ## x@MAFpopMask[popkey] <- value
                   change(cutoffs(x)$maxMAF, "popmask") <- value

                   x
                 })

setMethod("maxMAF", signature(x="VariantFilteringResults"),
          function(x) {
            warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
            if (!"MafDb" %in% param(x)$otherAnnotationsClass)
              stop("A MafDb object was not used to annotate variants.")

            cutoffs(x)[["maxMAF"]]$maxvalue
          })

setReplaceMethod("maxMAF", signature(x="VariantFilteringResults", value="numeric"),
                 function(x, value) {
                   warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
                   if (!"MafDb" %in% param(x)$otherAnnotationsClass)
                     stop("A MafDb object was not used to annotate variants.")

                   ## x@maxMAF <- value
                   change(cutoffs(x)$maxMAF, "maxvalue") <- value
                   x
                 })

setMethod("minPhastCons", signature(x="VariantFilteringResults"),
          function(x) {
            warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
            if (is.na(match("GScores", param(x)$otherAnnotationsClass)))
              stop("A GScores object was not used to annotate variants.")

            cutoffs(x)[["minPhastCons"]]$value
          })

setReplaceMethod("minPhastCons", signature(x="VariantFilteringResults", value="ANY"),
                 function(x, value) {
                   warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
                   if (is.na(match("GScores", param(x)$otherAnnotationsClass)))
                     stop("A phastConsDb object was not used to annotate variants.")

                   if (!is.na(value) && !is.numeric(value) && !is.integer(value))
                     stop("Only a numeric or NA value is allowed as minimum cutoff for phastCons scores.")

                   if (!is.na(value)) {
                     if (value < 0 || value > 1)
                       stop("The minimum cutoff for phastCons scores should be a number between 0 and 1.")
                   }

                   value <- as.numeric(value)
                   ## x@minPhastCons <- value
                   change(cutoffs(x)$phastCons100way, "value") <- value

                   x
                 })

setMethod("minPhylostratum", signature(x="VariantFilteringResults"),
          function(x) {
            warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
            if (is.na(match("GenePhylostrataDb", param(x)$otherAnnotationsClass)))
              stop("A GenePhylostrataDb object was not used to annotate variants.")

            ## x@minPhylostratumIndex
            cutoffs(x)[["genePhyloStratum"]]
          })

setReplaceMethod("minPhylostratum", signature(x="VariantFilteringResults", value="ANY"),
                 function(x, value) {
                   warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
                   whGPSdb <- match("GenePhylostrataDb", param(x)$otherAnnotationsClass)
                   if (is.na(whGPSdb))
                     stop("A GenePhylostrataDb object was not used to annotate variants.")

                   if (!is.na(value) && !is.numeric(value) && !is.integer(value) && !is.character(value))
                     stop("Only an integer, character or NA value is allowed as minimum phylostratum index.")

                   if (!is.na(value)) {
                     if (is.character(value)) {
                       value <- match(tolower(value),
                                      tolower(genePhylostrata(get(param(x)$otherAnnotations[whGPSdb]))$Description))
                       if (is.na(value))
                         stop(sprintf("%s is not valid. The minimum phylostratum character value should be one of: %s",
                                      value, paste(genePhylostrata(get(param(x)$otherAnnotations[whGPSdb]))$Description,
                                                   collapse=",")))
                     }

                     value <- as.integer(value)
                     if (value < 1 || value > nrow(genePhylostrata(get(param(x)$otherAnnotations[whGPSdb]))))
                       stop(sprintf("%d is not valid. The minimum phylostratum integer value should be between 1 and %d",
                                    nrow(genePhylostrata(get(param(x)$otherAnnotations[whGPSdb])))))
                   }

                   value <- as.integer(value)
                   ## x@minPhylostratumIndex <- value
                   change(cutoffs(x), "genePhyloStratum") <- value

                   x
                 })

setMethod("minScore5ss", signature(x="VariantFilteringResults"),
          function(x) {
            warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
            if (any(is.na(param(x)$spliceSiteMatricesFilenames)))
              stop("No splice site matrix was used to annotate variants.")

            x@minScore5ss
          })

setReplaceMethod("minScore5ss", signature(x="VariantFilteringResults", value="ANY"),
                 function(x, value) {
                   warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
                   if (any(is.na(param(x)$spliceSiteMatricesFilenames)))
                     stop("No splice site matrix was used to annotate variants.")

                   if (!is.na(value) && !is.numeric(value) && !is.integer(value))
                     stop("Only a numeric or NA value is allowed as minimum cutoff for cryptic 5'ss.")
                   x@minScore5ss <- as.numeric(value)
                   x
                 })

setMethod("minScore3ss", signature(x="VariantFilteringResults"),
          function(x) {
            warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
            if (any(is.na(param(x)$spliceSiteMatricesFilenames)))
              stop("No splice site matrix was used to annotate variants.")

            x@minScore3ss
          })

setReplaceMethod("minScore3ss", signature(x="VariantFilteringResults", value="ANY"),
                 function(x, value) {
                   warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
                   if (any(is.na(param(x)$spliceSiteMatricesFilenames)))
                     stop("No splice site matrix was used to annotate variants.")

                   if (!is.na(value) && !is.numeric(value) && !is.integer(value))
                     stop("Only a numeric or NA value is allowed as minimum cutoff for cryptic 3'ss.")
                   x@minScore3ss <- as.numeric(value)
                   x
                 })

setMethod("minCUFC", signature(x="VariantFilteringResults"),
          function(x) {
            warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
            x@minCUFC
          })

setReplaceMethod("minCUFC", signature(x="VariantFilteringResults", value="numeric"),
                 function(x, value) {
                   warning("This method is deprecated and has been replaced with the use of methods 'filters(x), active(x), cutoffs(x) and change(x). It will become defunct in the next release of Bioconductor.")
                   x@minCUFC <- value
                   x
                 })
