##
## show methodss
##

## provide a concise show method for CompressedVRangesList objects
## that result from doing split(vr, sampleNames(vr)) on a VRanges object 'vr'
setMethod("show", signature(object="CompressedVRangesList"),
          function(object) {
            lo <- length(object)
            cat(classNameForDisplay(object), " of length ", lo, "\n",
                sep = "")
            if (!is.null(names(object)))
              cat(S4Vectors:::labeledLine("names", names(object)))
          })

setMethod("show", signature(object="VariantFilteringResults"),
          function(object) {
            cat("\nVariantFiltering results object\n")
            cat(sprintf("\n  Genome version(s):"))
            for (i in seq_along(param(object)$seqInfos))
              if (length(param(object)$seqInfos[[i]]) > 0)
                cat(sprintf(" %s(%s)", paste(unique(genome(param(object)$seqInfos[[i]])), collapse=","),
                            seqlevelsStyle(param(object)$seqInfos[[i]])))

            sampleNames <- ifelse(length(samples(object)) <= 4, paste(samples(object), collapse=", "),
                                  paste(paste(head(samples(object), n=4), collapse=", "), "...", sep=", "))
            cat(sprintf("\n  Number of individuals: %d (%s)\n", length(samples(object)), sampleNames))

            if (inheritanceModel(object) != "any")
              cat(sprintf("  Variants segregate according to a(n) %s inheritance model\n", inheritanceModel(object)))
            else
              cat("  Variants are not filtered by inheritance model\n")
            if (length(param(object)@qualityFilterNames) > 0) {
              cat("  Quality filters\n")
              print(active(filters(object)[param(object)@qualityFilterNames]))
            }
            cat("  Functional annotation filters\n")
            print(active(filters(object))[setdiff(names(filters(object)), param(object)@qualityFilterNames)])
            if (minCUFC(object) > 0)
              cat(sprintf("    Minimum codon-usage abs log2-fold change: %s\n", minCUFC(object)))
            ## cat(sprintf("    Amino acid change type: %s\n", aaChangeType(object)))
            if (!all(variantLocation(object)))
              cat(sprintf("    Location restricted to: %s\n", paste(names(variantLocation(object))[variantLocation(object)], collapse=", ")))
            if (!all(variantConsequence(object)))
              cat(sprintf("    Consequence restricted to: %s\n", paste(names(variantConsequence(object))[variantConsequence(object)], collapse=", ")))
            if ("MafDb" %in% sapply(param(object)$otherAnnotations, class)) {
              cat(sprintf("    Populations used for MAF filtering: %s\n", paste(names(MAFpop(object))[MAFpop(object)], collapse=", ")))
              cat(sprintf("    Include MAF NA values: %s\n", ifelse(naMAF(object), "yes", "no")))
              cat(sprintf("    Maximum MAF: %.2f\n", maxMAF(object)))
            }
            if ("GScores" %in% sapply(param(object)$otherAnnotations, class)) {
              if (!is.na(minPhastCons(object)))
                cat(sprintf("    Minimum score for phastCons nucleotide conservation: %.2f\n", minPhastCons(object)))
            }
            if ("GenePhylostrataDb" %in% sapply(param(object)$otherAnnotations, class)) {
              whGPSdb <- match("GenePhylostrataDb", sapply(param(object)$otherAnnotations, class))
              if (!is.na(minPhylostratum(object)))
                cat(sprintf("    Minimum conserved gene phylostratum: %s (%d/%d)\n",
                            genePhylostrata(param(object)$otherAnnotations[[whGPSdb]])$Description[minPhylostratum(object)],
                            minPhylostratum(object),
                            nrow(genePhylostrata(param(object)$otherAnnotations[[whGPSdb]]))))
            }
            ## if (all(!is.na(param(object)$spliceSiteMatricesFilenames))) {
            ##   if (!is.na(minScore5ss(object)))
            ##     cat(sprintf("    Minimum score for cryptic 5'ss: %.2f\n", minScore5ss(object)))
            ##   if (!is.na(minScore3ss(object)))
            ##     cat(sprintf("    Minimum score for cryptic 3'ss: %.2f\n", minScore3ss(object)))
            ## }
          })

setMethod("summary", signature(object="VariantFilteringResults"),
          function(object, method=c("SO", "SOfull", "bioc")) {

            method <- match.arg(method)

            fv <- filteredVariants(object)
            fv1 <- fv[[1]]
            nvars <- length(unique(fv1$VCFIDX))

            res <- data.frame()
            if (method == "bioc") {
              total <- length(fv1)
              vxloc <- table(fv1$LOCATION)
              vxcon <- table(fv1$CONSEQUENCE)
              vxloc <- table(unique(mcols(fv1)[, c("VCFIDX", "LOCATION")])$LOCATION)
              vxcon <- table(unique(mcols(fv1)[, c("VCFIDX", "CONSEQUENCE")])$CONSEQUENCE)
              res <- data.frame("BIOCID"=c(names(vxloc), names(vxcon)),
                                "Nr. Variants"=c(as.integer(vxloc), as.integer(vxcon)),
                                check.names=FALSE)
              res <- res[res[[2]] > 0, ]
              f <- res[[2]] / nvars
              nd <- max(-floor(log10(f[f > 0]))-1)
              res$"% Variants" <- round(100*f, digits=nd)
            } else if (method == "SO" || method == "SOfull") {
              if (!"SOterms" %in% names(cutoffs(object)))
                stop("There is no 'SOterms' entry in the list of cutoff values.")

              soterms <- cutoffs(object)$SOterms
              gSO <- sog(object)

              ## no SO terms specified or no active SOterms filter, imply no restriction
              if (is.null(soterms) || all(is.na(soterms)) || !active(filters(object))["SOterms"])
                soterms <- nodes(gSO)[sapply(nodeData(gSO, nodes(gSO), "vcfIdx"), length) > 0]
              else
                soterms <- .findSOIDs(soterms, gSO)

              description <- character(0)
              numvariants <- integer(0)
              if (length(soterms) > 0) {
                if (method == "SO") { ## show only given SO terms
                  soterms <- soterms[order(match(soterms, tsort(gSO)))] ## show SO terms in topological order
                  numvariants <- sapply(nodeData(gSO, soterms, "vcfIdx"),
                                        function(vcfidx, allvcfidx) length(unique(vcfidx[vcfidx %in% allvcfidx])),
                                        fv1$VCFIDX)
                  description <- unlist(nodeData(gSO, soterms, "label"))
                  res <- data.frame("SOID"=soterms[numvariants > 0],
                                    "Description"=description[numvariants > 0],
                                    "Nr. Variants"=numvariants[numvariants > 0],
                                    check.names=FALSE, row.names=NULL, stringsAsFactors=FALSE)
                } else { ## show hierarchy involving the given SO terms

                  asoterms <- unique(c(soterms, .ancestorsSO(param(object), soterms)))
                  avcfidx <- nodeData(gSO, asoterms, "vcfIdx")
                  asoterms <- names(avcfidx)[sapply(avcfidx, length) > 0]
                  soterms <- unique(c(soterms, asoterms))
                  subgSO <- sog(param(object))
                  nodeDataDefaults(subgSO, "vcfIdx") <- integer(0)
                  nodeData(subgSO, soterms, "vcfIdx") <- nodeData(gSO, soterms, "vcfIdx")

                  dsoterms <- unique(c(soterms, .descendantsSO(param(object), soterms)))
                  subgSO <- subGraph(dsoterms, subgSO)
                  to <- tsort(subgSO)
                  soterms <- character(0)
                  for (v in to) {
                    soterms <- c(soterms, v)
                    parentterms <- inEdges(subgSO)[[v]]
                    parentvcfidx <- unlist(lapply(parentterms, function(x) nodeData(subgSO, x, "vcfIdx")), use.names=FALSE)
                    nodeData(subgSO, v, "vcfIdx")[[v]] <- unique(c(nodeData(subgSO, v, "vcfIdx")[[1]], parentvcfidx))
                    vcfidx <- nodeData(subgSO, v, "vcfIdx")[[v]]
                    numvariants <- c(numvariants, length(unique(vcfidx[vcfidx %in% fv1$VCFIDX])))
                    description <- c(description, unlist(nodeData(subgSO, v, "label"), use.names=FALSE))
                  }
                  alld <- dijkstra.sp(g=as(t(as(subgSO, "matrix")), "graphNEL"), start=to[length(to)])$distances
                  res <- data.frame("SOID"=soterms,
                                    "Level"=alld[soterms],
                                    "Description"=description,
                                    "Nr. Variants"=numvariants,
                                    check.names=FALSE, row.names=NULL, stringsAsFactors=FALSE)
                }
              }

              if (nrow(res) > 0) {
                f <- res[["Nr. Variants"]] / nvars
                nd <- max(-floor(log10(f[f > 0]))-1)
                res$"% Variants" <- round(100*f, digits=nd)
              } else
                res$"% Variants" <- numeric(0)
            }

            res
          })


##
## getter and setter methods
##

setMethod("length", "VariantFilteringResults", function(x) length(x@variants))

setMethod("filters", signature(x="VariantFilteringResults"),
          function(x) {
            x@filters
          })

setReplaceMethod("filters", signature(x="VariantFilteringResults"),
                 function(x, value) {
                   x@filters <- value
                   x
                 })

setMethod("cutoffs", signature(x="VariantFilteringResults"),
          function(x) {
            x@cutoffs
          })

setReplaceMethod("cutoffs", signature(x="VariantFilteringResults"),
                 function(x, value) {
                   x@cutoffs <- value
                   x
                 })

setMethod("softFilterMatrix", signature(x="VariantFilteringResults"),
          function(x) {
            softFilterMatrix(allVariants(x, groupBy="nothing"))
          })

setReplaceMethod("softFilterMatrix", signature(x="VariantFilteringResults"),
                 function(x, value) {
                   softFilterMatrix(x@variants) <- value
                   x
                 })

setMethod("sog", signature(x="VariantFilteringResults"),
          function(x) {
            x@gSO
          })

setMethod("samples", signature(object="VariantFilteringResults"),
          function(object) {
            object@activeSamples
          })

setReplaceMethod("samples", signature(object="VariantFilteringResults"),
                 function(object, value) {
                   if (inheritanceModel(object) != "unrelated individuals")
                     stop("Active samples can only be changed when no inheritance model is used (unrelated individuals).")

                   mask <- !value %in% param(object)$sampleNames
                   if (any(mask)) {
                     if (sum(mask) > 1)
                       stop(sprintf("%s are not valid sample names.", value[mask]))
                     else
                       stop(sprintf("%s is not a valid sample name.", value[mask]))
                   }

                   object@activeSamples <- value

                   object
                 })

setMethod("resetSamples", signature(object="VariantFilteringResults"),
          function(object) {
            object@activeSamples <- param(object)$sampleNames

            object
          })

setMethod("bamFiles", signature(object="VariantFilteringResults"),
          function(object) {
            object@bamViews
          })

setReplaceMethod("bamFiles", signature(object="VariantFilteringResults", value="BamViews"),
                 function(object, value) {
                   mask <- rownames(bamSamples(value)) %in% param(object)$sampleNames
                   if (!all(mask))
                     stop("Sample names do not match.")

                   object@bamViews <- value

                   object
                 })

setMethod("param", signature(x="VariantFilteringResults"),
          function(x) {
            x@inputParameters
          })

setMethod("inheritanceModel", signature(x="VariantFilteringResults"),
          function(x) {
            x@inheritanceModel
          })

setMethod("annoGroups", signature(x="VariantFilteringResults"),
          function(x) {
            x@annoGroups
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
            x@variantTypeMask
          })

setReplaceMethod("variantType", signature(x="VariantFilteringResults", value="logical"),
                 function(x, typkey=NA, value) {
                   if (any(is.na(value)))
                     stop("The given value(s) must be either TRUE or FALSE.")

                   if (is.na(typkey)) {
                     typkey <- names(x@variantTypeMask)
                     if (is.null(names(value))) {
                       if (length(value) > 1)
                         stop("When given multiple values they must be a logical vector whose names match the available location keywords.")
                       value <- do.call("names<-", list(rep(value, length(typkey)), typkey))
                     } else if (any(is.na(match(names(value), names(x@variantTypeMask)))))
                       stop(sprintf("Element names %s do not match the available variant type keywords",
                                    names(value)[is.na(match(names(value), names(x@variantTypeMask)))]))
                     else
                       typkey <- names(value)
                   } else {
                     if (any(is.na(match(typkey, names(x@variantTypeMask)))))
                       stop(sprintf("%s does not match the available variant type keywords", typkey))
                   }

                   x@variantTypeMask[typkey] <- value
                   cutoffs(x)$variantType <- x@variantTypeMask
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
                   cutoffs(x)$aaChangeType <- value
                   x
                 })

setMethod("variantLocation", signature(x="VariantFilteringResults"),
          function(x) {
            x@locationMask
          })

setReplaceMethod("variantLocation", signature(x="VariantFilteringResults", value="logical"),
                 function(x, lockey=NA, value) {
                   if (any(is.na(value)))
                     stop("The given value(s) must be either TRUE or FALSE.")

                   if (is.na(lockey)) {
                     lockey <- names(x@locationMask)
                     if (is.null(names(value))) {
                       if (length(value) > 1)
                         stop("When given multiple values they must be a logical vector whose names match the available location keywords.")
                       value <- do.call("names<-", list(rep(value, length(lockey)), lockey))
                     } else if (any(is.na(match(names(value), names(x@locationMask)))))
                       stop(sprintf("Element names %s do not match the available location keywords",
                                    names(value)[is.na(match(names(value), names(x@locationMask)))]))
                     else
                       lockey <- names(value)
                   } else {
                     if (any(is.na(match(lockey, names(x@locationMask)))))
                       stop(sprintf("%s does not match the available location keywords", lockey))
                   }

                   x@locationMask[lockey] <- value
                   x
                 })

setMethod("variantConsequence", signature(x="VariantFilteringResults"),
          function(x) {
            x@consequenceMask
          })

setReplaceMethod("variantConsequence", signature(x="VariantFilteringResults", value="logical"),
                 function(x, conkey=NA, value) {

                   if (any(is.na(value)))
                     stop("The given value(s) must be either TRUE or FALSE.")

                   if (is.na(conkey)) {
                     conkey <- names(x@consequenceMask)
                     if (is.null(names(value))) {
                       if (length(value) > 1)
                         stop("When given multiple values they must be a logical vector whose names match the available consequence keywords.")
                       value <- do.call("names<-", list(rep(value, length(conkey)), conkey))
                     } else if (any(is.na(match(names(value), names(x@consequenceMask)))))
                       stop(sprintf("Element names %s do not match the available consequence keywords",
                                    names(value)[is.na(match(names(value), names(x@consequenceMask)))]))
                     else
                       conkey <- names(value)
                   } else {
                     if (any(is.na(match(conkey, names(x@consequenceMask)))))
                       stop(sprintf("%s does not match the available consequence keywords", conkey))
                   }

                   x@consequenceMask[conkey] <- value
                   x
                 })

setMethod("naMAF", signature(x="VariantFilteringResults"),
          function(x) {
            if (!any(c("MafDb", "MafDb2") %in% sapply(param(x)$otherAnnotations, class)))
              stop("A MafDb object was not used to annotate variants.")

            x@naMAF
          })

setReplaceMethod("naMAF", signature(x="VariantFilteringResults", value="logical"),
          function(x, value) {
            if (!any(c("MafDb", "MafDb2") %in% sapply(param(x)$otherAnnotations, class)))
              stop("A MafDb object was not used to annotate variants.")

            x@naMAF <- value
            x
          })

setMethod("MAFpop", signature(x="VariantFilteringResults"),
          function(x) {
            if (!any(c("MafDb", "MafDb2") %in% sapply(param(x)$otherAnnotations, class)))
              stop("A MafDb object was not used to annotate variants.")

            x@MAFpopMask
          })

setReplaceMethod("MAFpop", signature(x="VariantFilteringResults", value="logical"),
                 function(x, popkey=NA, value) {
                   if (!any(c("MafDb", "MafDb2") %in% sapply(param(x)$otherAnnotations, class)))
                     stop("A MafDb object was not used to annotate variants.")

                   if (any(is.na(value)))
                     stop("The given value(s) must be either TRUE or FALSE.")

                   if (is.na(popkey)) {
                     popkey <- names(x@MAFpopMask)
                     if (is.null(names(value))) {
                       if (length(value) > 1)
                         stop("When given multiple values they must be a logical vector whose names match the available population keywords.")
                       value <- do.call("names<-", list(rep(value, length(popkey)), popkey))
                     } else if (any(is.na(match(names(value), names(x@MAFpopMask)))))
                       stop(sprintf("Element names %s do not match the available population keywords",
                                    names(value)[is.na(match(names(value), names(x@MAFpopMask)))]))
                     else
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
            if (!any(c("MafDb", "MafDb2") %in% sapply(param(x)$otherAnnotations, class)))
              stop("A MafDb object was not used to annotate variants.")

            x@maxMAF
          })

setReplaceMethod("maxMAF", signature(x="VariantFilteringResults", value="numeric"),
                 function(x, value) {
                   if (!any(c("MafDb", "MafDb2") %in% sapply(param(x)$otherAnnotations, class)))
                     stop("A MafDb object was not used to annotate variants.")

                   x@maxMAF <- value
                   x
                 })

setMethod("minPhastCons", signature(x="VariantFilteringResults"),
          function(x) {
            if (is.na(match("GScores", sapply(param(x)$otherAnnotations, class))))
              stop("A GScores object was not used to annotate variants.")

            x@minPhastCons
          })

setReplaceMethod("minPhastCons", signature(x="VariantFilteringResults", value="ANY"),
                 function(x, value) {
                   if (is.na(match("GScores", sapply(param(x)$otherAnnotations, class))))
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

setMethod("minScore5ss", signature(x="VariantFilteringResults"),
          function(x) {
            if (any(is.na(param(x)$spliceSiteMatricesFilenames)))
              stop("No splice site matrix was used to annotate variants.")

            x@minScore5ss
          })

setReplaceMethod("minScore5ss", signature(x="VariantFilteringResults", value="ANY"),
                 function(x, value) {
                   if (any(is.na(param(x)$spliceSiteMatricesFilenames)))
                     stop("No splice site matrix was used to annotate variants.")

                   if (!is.na(value) && !is.numeric(value) && !is.integer(value))
                     stop("Only a numeric or NA value is allowed as minimum cutoff for cryptic 5'ss.")
                   x@minScore5ss <- as.numeric(value)
                   x
                 })

setMethod("minScore3ss", signature(x="VariantFilteringResults"),
          function(x) {
            if (any(is.na(param(x)$spliceSiteMatricesFilenames)))
              stop("No splice site matrix was used to annotate variants.")

            x@minScore3ss
          })

setReplaceMethod("minScore3ss", signature(x="VariantFilteringResults", value="ANY"),
                 function(x, value) {
                   if (any(is.na(param(x)$spliceSiteMatricesFilenames)))
                     stop("No splice site matrix was used to annotate variants.")

                   if (!is.na(value) && !is.numeric(value) && !is.integer(value))
                     stop("Only a numeric or NA value is allowed as minimum cutoff for cryptic 3'ss.")
                   x@minScore3ss <- as.numeric(value)
                   x
                 })

setMethod("minCUFC", signature(x="VariantFilteringResults"),
          function(x) {
            x@minCUFC
          })

setReplaceMethod("minCUFC", signature(x="VariantFilteringResults", value="numeric"),
                 function(x, value) {
                   x@minCUFC <- value
                   x
                 })


setMethod("filters", signature(x="VariantFilteringResults"),
          function(x) {
            x@filters
          })


##
## methods for retrieving variants
##

## get all variants without applying any filter
setMethod("allVariants", signature(x="VariantFilteringResults"), 
          function(x, groupBy="sample") {
            vars <- x@variants
            if (groupBy[1] %in% "sample") {
              f <- sampleNames(x@variants)
              if (all(is.na(f)))
                f <- rep("nosample", length(x))
              vars <- split(x@variants, f)
            } else if (groupBy[1] %in% colnames(mcols(x@variants)))
              vars <- split(x@variants, mcols(x@variants)[, groupBy])

            vars
          })

## get variants after applying all filters
setMethod("filteredVariants", signature(x="VariantFilteringResults"), 
          function(x, groupBy="sample", unusedColumns.rm=FALSE) {
            if (length(x@filters) > 0)
              x <- softFilter(x, filters(x))
            varsxsam <- allVariants(x)
            vars <- varsxsam[[1]]
            selcols <- names(active(filters(x)))[active(filters(x))] ## default VCF filters may or may not show up
            rowsMask <- apply(softFilterMatrix(vars)[, selcols, drop=FALSE], 1, all, na.rm=TRUE)
            if (!all(param(x)$sampleNames %in% samples(x)) && all(samples(x) %in% names(varsxsam))) { ## not all samples are active
              varsxsam <- varsxsam[samples(x)]

              ## discard variants that are not present in active samples
              gt <- sapply(varsxsam, function(x) x$GT)
              gt <- apply(gt, 1, paste, collapse="")
              gt <- gsub("/|\\|", "", gt)
              gt <- gsub(".", "0", gt, fixed=TRUE)
              gt <- gsub("0+", "0", gt)
              rowsMask <- rowsMask & !gt %in% "0"
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
            if (!all(variantType(x)))
              rowsMask <- rowsMask & vars$TYPE %in% names(variantType(x))[variantType(x)]

            ## location of variant
            if (!all(variantLocation(x)))
              rowsMask <- rowsMask & vars$LOCATION %in% names(variantLocation(x))[variantLocation(x)]

            ## consquence of variant
            if (!all(variantConsequence(x)))
              rowsMask <- rowsMask & vars$CONSEQUENCE %in% names(variantConsequence(x))[variantConsequence(x)]

            ## type of amino acid change
            if (aaChangeType(x) != "Any")
              rowsMask <- rowsMask & vars$AAchangeType == aaChangeType(x)

            ## minimum allele frequency
            mtNoMAF <- NULL
            if ("MafDb" %in% sapply(param(x)$otherAnnotations, class)) {
              mtNoMAF <- match(names(MAFpop(x)), colnames(mcols(vars)))
              maxMAFannot <- rep(NA_real_, length(vars))
              if (any(MAFpop(x))) {
                mtNoMAF <- NULL
                maxMAFannot <- do.call(pmax, c(as.list(mcols(vars[, names(MAFpop(x))[MAFpop(x)]])), na.rm=TRUE))
                if (naMAF(x))
                  maxMAFannot[is.na(maxMAFannot)] <- -Inf
                else
                  maxMAFannot[is.na(maxMAFannot)] <- Inf

                rowsMask <- rowsMask & maxMAFannot <= maxMAF(x)
                rowsMask[is.na(rowsMask)] <- FALSE

                maxMAFannot[!is.finite(maxMAFannot)] <- NA_real_
                if (any(!MAFpop(x)))
                  mtNoMAF <- match(names(MAFpop(x))[!MAFpop(x)], colnames(mcols(vars)))
              }
            }

            ## nucleotide conservation
            mtNoMinPhastCons <- NULL
            if (!is.na(match("GScores", sapply(param(x)$otherAnnotations, class)))) {
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
            ## mtNoSCORE5ss <- mtNoSCORE3ss <- NULL
            ## if (all(!is.na(param(x)$spliceSiteMatricesFilenames))) {
            ##   crypssMask <- rep(FALSE, length(vars))
            ##   if (is.na(minScore5ss(x)))
            ##     mtNoSCORE5ss <- grep("SCORE5ss", colnames(mcols(vars)))
            ##   else {
            ##     cryp5ssMask <- vars$SCORE5ssALT >= minScore5ss(x)
            ##     cryp5ssMask[is.na(cryp5ssMask)] <- FALSE
            ##     crypssMask <- crypssMask | cryp5ssMask
            ##   }
            ##   if (is.na(minScore3ss(x)))
            ##     mtNoSCORE3ss <- grep("SCORE3ss", colnames(mcols(vars)))
            ##   else {
            ##     cryp3ssMask <- vars$SCORE3ssALT >= minScore3ss(x)
            ##     cryp3ssMask[is.na(cryp3ssMask)] <- FALSE
            ##     crypssMask <- crypssMask | cryp3ssMask
            ##   }
            ##   ## if no filter on 5' and 3' cryptic ss is set, then select all rows
            ##   if (is.na(minScore5ss(x)) && is.na(minScore3ss(x)))
            ##     crypssMask <- rep(TRUE, length(vars))
            ##
            ##   rowsMask <- rowsMask & crypssMask
            ## }

            ## codon-usage fold-change
            minCUFCannot <- log2(vars$CUALT) - log2(vars$CUREF)
            naCUFCmask <- rep(TRUE, length(vars))
            minCUFCannot[is.na(minCUFCannot)] <- Inf

            rowsMask <- rowsMask & abs(minCUFCannot) >= minCUFC(x)
            rowsMask[is.na(rowsMask)] <- FALSE

            minCUFCannot[!is.finite(minCUFCannot)] <- NA_real_

            ## select variant annotations using the logical mask of the filters
            colsIdx <- setdiff(1:ncol(mcols(vars)), mtNoMAF)
            if (unusedColumns.rm) ## remove data columns that are not used for filtering
              colsIdx <- setdiff(colsIdx, c(mtNoMinPhastCons, mtNoMinPhylostratum))## , mtNoSCORE5ss, mtNoSCORE3ss))

            colsMask <- rep(FALSE, ncol(mcols(vars)))
            colsMask[colsIdx] <- TRUE

            vars <- unlist(varsxsam, use.names=FALSE)
            sampleNames(vars) <- Rle(names(varsxsam), elementNROWS(varsxsam))
            rowsMask <- rep(rowsMask, length(varsxsam))
            if ("MafDb" %in% sapply(param(x)$otherAnnotations, class))
              vars$maxMAF <- rep(maxMAFannot, length(varsxsam))
            vars$CUFC <- rep(minCUFCannot, length(varsxsam))

            vars <- vars[rowsMask, colsMask]

            if (groupBy[1] %in% "sample")
              vars <- split(vars, sampleNames(vars))
            else if (groupBy[1] %in% colnames(mcols(vars)))
              vars <- split(vars, mcols(vars)[, groupBy])

            vars
          })

##
## shiny app to filter and visualize variants
##

setMethod("reportVariants", signature(vfResultsObj="VariantFilteringResults"),
          function(vfResultsObj, type=c("shiny", "csv", "tsv"), file=NULL, UCSCorg=NA_character_) {
              
  type <- match.arg(type)

  if (class(vfResultsObj) != "VariantFilteringResults")
    stop("Input argument 'vfResultsObj' should be an object of class 'VariantFilteringResults'.")

  if ((type == "csv" || type == "tsv") && is.null(file))
    stop("If type=\"csv\" or type=\"tsv\" then the input argument 'file' cannot be NULL.")

  if (type == "csv" || type == "tsv") {
    varsdf <- filteredVariants(vfResultsObj, groupBy="nothing")
    varsdf <- as.data.frame(DataFrame(CHR=seqnames(varsdf),
                                      POS=start(varsdf),
                                      SAMPLEID=sampleNames(varsdf),
                                      mcols(varsdf)))
    firstcols <- c("VARID", "dbSNP", "CHR", "POS", "SAMPLEID", "GT")
    varsdf <- varsdf[, c(firstcols, setdiff(colnames(varsdf), firstcols))]
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
  ## the 'organism()' getter for BSgenome objects when not specified with the
  ## function argument 'UCSCorg'

  if (is.na(UCSCorg)) {
    UCSCorg <- organism(param(vfResultsObj)$bsgenome)
    UCSCorg <- paste0(substr(UCSCorg, 1, 1), strsplit(UCSCorg, " ")[[1]][2])
  }
  genomeBuild <- providerVersion(param(vfResultsObj)$bsgenome)

  annotationObjClasses <- sapply(param(vfResultsObj)$otherAnnotations, class)
  mtMafDb <- match("MafDb", annotationObjClasses)
  mtPhastConsDb <- match("GScores", annotationObjClasses)
  mtGenePhylostrataDb <- match("GenePhylostrataDb", annotationObjClasses)
  phylostrata <- "Unavailable"
  if (!is.na(mtGenePhylostrataDb))
    phylostrata <- rev(genePhylostrata(param(vfResultsObj)$otherAnnotations[[mtGenePhylostrataDb]])$Description)
  ## crypsplice <- param(vfResultsObj)$spliceSiteMatricesFilenames
  
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
                conditionalPanel(condition="input.tsp == 'genome'", helpText(strong("Restrict variants to the following types:"))),
                conditionalPanel(condition="input.tsp == 'genome'", checkboxInput("SNV", "SNV", TRUE)),
                conditionalPanel(condition="input.tsp == 'genome'", checkboxInput("Insertion", "Insertion", TRUE)),
                conditionalPanel(condition="input.tsp == 'genome'", checkboxInput("Deletion", "Deletion", TRUE)),
                conditionalPanel(condition="input.tsp == 'genome'", checkboxInput("MNV", "MNV", TRUE)),
                conditionalPanel(condition="input.tsp == 'genome'", checkboxInput("Delins", "Delins", TRUE)),
                ## transcript tab
                conditionalPanel(condition="input.tsp == 'transcript'", numericInput('minCUFC', 'Minimum Codon Usage Absolute log2-Fold Change:', 0.00)),
                conditionalPanel(condition="input.tsp == 'transcript'", helpText(strong("Restrict variants to the following locations:"))),
                conditionalPanel(condition="input.tsp == 'transcript'", checkboxInput("coding", "Coding", TRUE)),
                conditionalPanel(condition="input.tsp == 'transcript'", checkboxInput("fiveUTR", "5' UTR", TRUE)),
                conditionalPanel(condition="input.tsp == 'transcript'", checkboxInput("threeUTR", "3' UTR", TRUE)),
                conditionalPanel(condition="input.tsp == 'transcript'", checkboxInput("intron", "Intronic", TRUE)),
                conditionalPanel(condition="input.tsp == 'transcript'", checkboxInput("spliceSite", "Known splice site", TRUE)),
                conditionalPanel(condition="input.tsp == 'transcript'", checkboxInput("promoter", "Promoter", TRUE)),
                conditionalPanel(condition="input.tsp == 'transcript'", checkboxInput("intergenic", "Intergenic", TRUE)),
                
                ## gene tab
                conditionalPanel(condition="input.tsp == 'gene'",
                                 checkboxInput('OMIMpresentFlag',
                                               'Filter by presence in OMIM', FALSE)),
                conditionalPanel(condition="input.tsp == 'gene' && input.OMIMpresentFlag == true",
                                 selectInput("OMIMpresent", "Present in OMIM:",
                                             choices=c("Yes", "No"))),
                ## protein tab
                conditionalPanel(condition="input.tsp == 'protein'", helpText(strong("Restrict variants to the following consequences:"))),
                conditionalPanel(condition="input.tsp == 'protein'", checkboxInput("synonymous", "Synonymous", TRUE)),
                conditionalPanel(condition="input.tsp == 'protein'", checkboxInput("nonsynonymous", "Nonsynonymous", TRUE)),
                conditionalPanel(condition="input.tsp == 'protein'", checkboxInput("frameshift", "Frameshift", TRUE)),
                conditionalPanel(condition="input.tsp == 'protein'", checkboxInput("nonsense", "Nonsense", TRUE)),
                conditionalPanel(condition="input.tsp == 'protein'", checkboxInput("not tranlated", "Not translated", TRUE)),
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
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('AFExAC', 'All MAF ExAC', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('AFR_AFExAC', 'AFR MAF ExAC', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('AMR_AFExAC', 'AMR MAF ExAC', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('Adj_AFExAC', 'Adj MAF ExAC', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('EAS_AFExAC', 'EAS MAF ExAC', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('FIN_AFExAC', 'FIN MAF ExAC', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('NFE_AFExAC', 'NFE MAF ExAC', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('OTH_AFExAC', 'OTH MAF ExAC', TRUE)),
                conditionalPanel(condition="input.tsp == 'maf'", checkboxInput('SAS_AFExAC', 'SAS MAF ExAC', TRUE)),
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
                conditionalPanel(condition="input.tsp == 'cryp'", checkboxInput('minScore5ssFlag', 'Filter by cryptic 5\'ss', FALSE)),
                conditionalPanel(condition="input.tsp == 'cryp' && input.minScore5ssFlag == true",
                                 numericInput('minScore5ss', 'Minimum cryptic 5\'ss score:', -99)),
                conditionalPanel(condition="input.tsp == 'cryp'", checkboxInput('minScore3ssFlag', 'Filter by cryptic 3\'ss', FALSE)),
                conditionalPanel(condition="input.tsp == 'cryp' && input.minScore3ssFlag == true",
                                 numericInput('minScore3ss', 'Minimum cryptic 3\'ss score:', -99)),
                tags$hr(),
                downloadButton('downloadData', 'Download Variants'),
                downloadButton('generateReport', 'Generate Report'),
                actionButton('closesavebutton', 'Save & Close')
              ),
              mainPanel(
                uiOutput('mytabs')
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
        varTypMask <- variantType(vfResultsObj)
        for (i in names(varTypMask))
          varTypMask[i] <- input[[i]]
        variantType(vfResultsObj) <- varTypMask

        ## variant location
        locMask <- variantLocation(vfResultsObj)
        for (i in names(locMask))
          locMask[i] <- input[[i]]
        variantLocation(vfResultsObj) <- locMask

        ## variant consequence
        conMask <- variantConsequence(vfResultsObj)
        for (i in names(conMask))
          conMask[i] <- input[[i]]
        variantConsequence(vfResultsObj) <- conMask

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
        ## if (all(!is.na(crypsplice))) {
        ##   if (input$minScore5ssFlag)
        ##     minScore5ss(vfResultsObj) <- input$minScore5ss
        ##   else
        ##     minScore5ss(vfResultsObj) <- NA_real_
        ##
        ##   if (input$minScore3ssFlag)
        ##     minScore3ss(vfResultsObj) <- input$minScore3ss
        ##   else
        ##     minScore3ss(vfResultsObj) <- NA_real_
        ## }

        ## codon-usage fold-change
        minCUFC(vfResultsObj) <- as.numeric(input$minCUFC)

        stopApp(returnValue=vfResultsObj)
      })
    })

    filteredVariantsReact <- reactive({
      vdf <- data.frame()
      withProgress(message="Filtering variants", value=0, {
      incProgress(1/3, detail="setting filters ...")
      ## presence in dbSNP
      ## if (input$dbSNPpresentFlag)
      ##   dbSNPpresent(vfResultsObj) <- input$dbSNPpresent
      ## else
      ##   dbSNPpresent(vfResultsObj) <- NA_character_

      ## presence in OMIM
      ## if (input$OMIMpresentFlag)
      ##   OMIMpresent(vfResultsObj) <- input$OMIMpresent
      ## else
      ##   OMIMpresent(vfResultsObj) <- NA_character_

      ## type of variant
      ## varTypMask <- variantType(vfResultsObj)
      ## for (i in names(varTypMask))
      ##   varTypMask[i] <- input[[i]]
      ## variantType(vfResultsObj) <- varTypMask

      ## variant location
      ## locMask <- variantLocation(vfResultsObj)
      ## for (i in names(locMask))
      ##   locMask[i] <- input[[i]]
      ## variantLocation(vfResultsObj) <- locMask

      ## variant consequence
      ## conMask <- variantConsequence(vfResultsObj)
      ## for (i in names(conMask))
      ##   conMask[i] <- input[[i]]
      ## variantConsequence(vfResultsObj) <- conMask

      ## type of amino acid change
      ## aaChangeType(vfResultsObj) <- input$aaChangeType

      ## minimum allele frequency
      ## if (!is.na(mtMafDb)) {
      ##   mafMask <- MAFpop(vfResultsObj)
      ##   for (i in names(mafMask)) ## unlisting a reactivevalues object does not work :(
      ##     mafMask[i] <- input[[i]]
      ##   MAFpop(vfResultsObj) <- mafMask
      ##   maxMAF(vfResultsObj) <- as.numeric(input$maxMAF)
      ##   naMAF(vfResultsObj) <- input$naMAF
      ## }

      ## nucleotide conservation
      ## if (!is.na(mtPhastConsDb)) {
      ##   if (input$minPhastConsFlag)
      ##     minPhastCons(vfResultsObj) <- input$minPhastCons
      ##   else
      ##     minPhastCons(vfResultsObj) <- NA_real_
      ## }

      ## gene conservation
      ## if (!is.na(mtGenePhylostrataDb)) {
      ##   if (input$minPhylostratumFlag)
      ##     minPhylostratum(vfResultsObj) <- input$minPhylostratum
      ##   else
      ##     minPhylostratum(vfResultsObj) <- NA_integer_
      ## }

      ## cryptic splice sites
      ## if (all(!is.na(crypsplice))) {
      ##   if (input$minScore5ssFlag)
      ##     minScore5ss(vfResultsObj) <- as.numeric(input$minScore5ss)
      ##   else
      ##     minScore5ss(vfResultsObj) <- NA_real_
      ##
      ##   if (input$minScore3ssFlag)
      ##     minScore3ss(vfResultsObj) <- as.numeric(input$minScore3ss)
      ##   else
      ##     minScore3ss(vfResultsObj) <- NA_real_
      ## }

      ## codon-usage fold-change
      ## minCUFC(vfResultsObj) <- as.numeric(input$minCUFC)

      incProgress(2/3, detail="applying filters ...")
      fvxsam <- filteredVariants(vfResultsObj)
      incProgress(3/3, detail="retrieving variants ...")
      fv <- fvxsam[[1]]
      vdf <- DataFrame(VarID=fv$VARID, CHR=seqnames(fv),
                       POS=start(fv), POSITION=start(fv),
                       mcols(fv))
      vdf$REF <- ref(fv)
      vdf$ALT <- alt(fv)
      if (length(fv) > 1) {
        vdf$DP <- as.integer(round(rowMeans(sapply(fvxsam, totalDepth)), digits=0))
        vdf$REFDP <- as.integer(round(rowMeans(sapply(fvxsam, refDepth)), digits=0))
        vdf$ALTDP <- as.integer(round(rowMeans(sapply(fvxsam, altDepth)), digits=0))
      } else if (length(fv) == 1) {
        vdf$DP <- as.integer(round(mean(sapply(fvxsam, totalDepth)), digits=0))
        vdf$REFDP <- as.integer(round(mean(sapply(fvxsam, refDepth)), digits=0))
        vdf$ALTDP <- as.integer(round(mean(sapply(fvxsam, altDepth)), digits=0))
      } else
        vdf$DP <- vdf$REFDP <- vdf$ALTDP <- integer()
      vdf <- as.data.frame(vdf)
      
      if (nrow(vdf) > 0) {
        varlocs <- paste0("<a href=http://genome.ucsc.edu/cgi-bin/hgTracks?org=", UCSCorg,
                          "&db=", genomeBuild, 
                          "&position=", vdf$CHR, ":", vdf$POS, " target=\"ucsc\">",
                          vdf$CHR, ":", vdf$POS,
                          "</a>")
        tempOMIM <- rep(NA_character_, nrow(vdf))
        tempOMIM[!is.na(vdf$OMIM)] <- sapply(strsplit(vdf$OMIM[!is.na(vdf$OMIM)], " *, *"), function(x) {
          z <- paste0("<a href=http://www.omim.org/entry/", x, " target=\"omim\">", x, "</a>")
          a <- paste(z, collapse=", ")
          a
        })
      
        tempdbSNP <- rep(NA_character_, nrow(vdf))
        tempdbSNP[!is.na(vdf$dbSNP)] <- sapply(strsplit(gsub("rs", "", vdf$dbSNP[!is.na(vdf$dbSNP)]), " *, *"),
                                               function(x) {
          z <- paste0("<a href=http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=", x,
                      " target=\"dbsnp\">", paste0("rs", x), "</a>")
          a <- paste(z, collapse=", ")
          a
        })
      
        vdf[["POSITION"]] <- varlocs
        vdf[["OMIM"]] <- tempOMIM
        vdf[["dbSNP"]] <- tempdbSNP
      }

      vdf <- vdf[, -match(c("CHR", "POS"), colnames(vdf))]
      rownames(vdf) <- NULL

      }) ## withProgress
      vdf 
    })

    output$mytabs <- renderUI({
      tabPanelList <- list(tabPanel("Genome", tableOutput('tableGenome'), value="genome"),
                           tabPanel("Gene", tableOutput('tableGene'), value="gene"),
                           tabPanel("Transcript", tableOutput('tableTranscript'), value="transcript"),
                           tabPanel("Protein", htmlOutput('tableProtein'), value="protein"))

      if ("MafDb" %in% annotationObjClasses)
        tabPanelList[[length(tabPanelList)+1]] <- tabPanel("MAF", tableOutput('tableMAF'), value="maf")

      if ("GScores" %in% annotationObjClasses || "GenePhylostrataDb" %in% annotationObjClasses)
        tabPanelList[[length(tabPanelList)+1]] <- tabPanel("Conservation", tableOutput('tableConservation'), value="conservation")

      ## if (all(!is.na(param(vfResultsObj)$spliceSiteMatricesFilenames)))
      ##   tabPanelList[[length(tabPanelList)+1]] <- tabPanel("CrypSplice", tableOutput('tableCrypSplice'), value="cryp")

      tabPanelList[[length(tabPanelList)+1]] <- tabPanelAbout()

      tabPanelList[[length(tabPanelList)+1]] <- "tsp"

      names(tabPanelList) <- c(rep("", length(tabPanelList)-1), "id")

      do.call(tabsetPanel, tabPanelList)
    })

    output$tableGenome <- renderTable({
      filteredVariantsReact()[, c("VarID", "POSITION", "dbSNP", "TYPE", "HGVSg", "DP", "REFDP", "ALTDP")]
    }, NA.string="NA",  sanitize.text.function=function(x){x})

    output$tableGene <- renderTable({
      filteredVariantsReact()[, c("VarID", "POSITION", "GENE", "LOCATION", "HGVSc", "OMIM")]
    }, NA.string="NA",  sanitize.text.function=function(x){x})

    output$tableTranscript <- renderTable({
      filteredVariantsReact()[, c("VarID", "POSITION", "GENE", "TXID", "LOCATION", "LOCSTART", "LOCEND", "cDNALOC", "HGVSc", "CUREF", "CUALT", "CUFC")]
    }, NA.string="NA",  sanitize.text.function=function(x){x})

    output$tableProtein <- renderTable({
      selcols <- c("VarID", "GENE", "HGVSp", "CONSEQUENCE", "AAchange", "AAchangeType")
      if ("PolyPhenDb" %in% annotationObjClasses)
        selcols <- c(selcols, "PolyPhen2")
      if ("PROVEANDb" %in% annotationObjClasses)
        selcols <- c(selcols, "PROVEAN")

      filteredVariantsReact()[, selcols]
    }, NA.string="NA",  sanitize.text.function=function(x){x})

    output$tableMAF <- renderTable({
      selcols <- selcolsnames <- c("VarID", "dbSNP")
      if ("MafDb" %in% annotationObjClasses) {
        selcols <- c(selcols, "maxMAF", names(MAFpop(vfResultsObj))[MAFpop(vfResultsObj)])
        selcolsnames <- c(selcolsnames, "Max", gsub("_", " ", names(MAFpop(vfResultsObj))[MAFpop(vfResultsObj)]))
      }

      fv <- filteredVariantsReact()
      fv <- fv[, selcols]
      colnames(fv) <- selcolsnames

      fv
    }, NA.string="NA",  sanitize.text.function=function(x){x})

    output$tableConservation <- renderTable({
      selcols <- c("VarID", "POSITION", "GENE")
      if ("GScores" %in% annotationObjClasses && !is.na(mtPhastConsDb)) {
        cname <- type(param(vfResultsObj)$otherAnnotations[[mtPhastConsDb]])
        selcols <- c(selcols, cname)
      }
      if ("GenePhylostrataDb" %in% annotationObjClasses)
        selcols <- c(selcols, "GenePhylostratum", "GenePhylostratumTaxID")

      fv <- filteredVariantsReact()[, selcols]

      if ("GenePhylostrataDb" %in% annotationObjClasses) {
        fv$GenePhylostratum[!is.na(fv$GenePhylostratum)] <- sprintf("<a href=\"http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s&lvl=5\" target='ncbitaxonomybrowser'>%s</a>", fv$GenePhylostratumTaxID[!is.na(fv$GenePhylostratum)], fv$GenePhylostratum[!is.na(fv$GenePhylostratum)])
        fv <- fv[, -match("GenePhylostratumTaxID", colnames(fv))]
      }

      fv
    }, NA.string="NA", sanitize.text.function=function(x){x})

    output$tableCrypSplice <- renderTable({
      selcols <- selcolsnames <- c("VarID", "POSITION")

      ## if (all(!is.na(param(vfResultsObj)$spliceSiteMatricesFilenames))) {
      ##   selcols <- c(selcols, "SCORE5ssREF", "SCORE5ssALT", "SCORE5ssPOS", "SCORE3ssREF", "SCORE3ssALT", "SCORE3ssPOS")
      ##   selcolsnames <- c(selcolsnames, "5'ss Ref", "5'ss Alt", "5'ss Pos", "3'ss Ref", "3'ss Alt", "3'ss Pos")
      ## }

      fv <- filteredVariantsReact()[, selcols]
      colnames(fv) <- selcolsnames

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
        varTypMask <- variantType(vfResultsObj)
        for (i in names(varTypMask))
          varTypMask[i] <- input[[i]]
        variantType(vfResultsObj) <- varTypMask
 
        ## variant location
        locMask <- variantLocation(vfResultsObj)
        for (i in names(locMask))
          locMask[i] <- input[[i]]
        variantLocation(vfResultsObj) <- locMask

        ## variant consequence
        conMask <- variantConsequence(vfResultsObj)
        for (i in names(conMask))
          conMask[i] <- input[[i]]
        variantConsequence(vfResultsObj) <- conMask

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
        ## if (all(!is.na(crypsplice))) {
        ##   if (input$minScore5ssFlag)
        ##     minScore5ss(vfResultsObj) <- as.numeric(input$minScore5ss)
        ##   else
        ##     minScore5ss(vfResultsObj) <- NA_real_
        ##
        ##   if (input$minScore3ssFlag)
        ##     minScore3ss(vfResultsObj) <- as.numeric(input$minScore3ss)
        ##   else
        ##     minScore3ss(vfResultsObj) <- NA_real_
        ## }

        ## codon-usage fold-change
        minCUFC(vfResultsObj) <- as.numeric(input$minCUFC)

        ## apply filters
        fvxsam <- filteredVariants(vfResultsObj)
        fv <- fvxsam[[1]]

        ## build data frame with the annotated and filtered variants
        vdf <- DataFrame(VarID=fv$VARID, CHR=seqnames(fv),
                         POS=start(fv), POSITION=start(fv),
                         mcols(fv))
        vdf$REF <- ref(fv)
        vdf$ALT <- alt(fv)
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
        varTypMask <- variantType(vfResultsObj)
        for (i in names(varTypMask))
          varTypMask[i] <- input[[i]]
        variantType(vfResultsObj) <- varTypMask

        ## variant location
        locMask <- variantLocation(vfResultsObj)
        for (i in names(locMask))
          locMask[i] <- input[[i]]
        variantLocation(vfResultsObj) <- locMask

        ## variant consequence
        conMask <- variantConsequence(vfResultsObj)
        for (i in names(conMask))
          conMask[i] <- input[[i]]
        variantConsequence(vfResultsObj) <- conMask

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
        ## if (all(!is.na(crypsplice))) {
        ##   if (input$minScore5ssFlag)
        ##     minScore5ss(vfResultsObj) <- as.numeric(input$minScore5ss)
        ##   else
        ##     minScore5ss(vfResultsObj) <- NA_real_
        ##
        ##   if (input$minScore3ssFlag)
        ##     minScore3ss(vfResultsObj) <- as.numeric(input$minScore3ss)
        ##   else
        ##     minScore3ss(vfResultsObj) <- NA_real_
        ## }

        ## codon-usage fold-change
        minCUFC(vfResultsObj) <- as.numeric(input$minCUFC)

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
        if (!all(variantType(vfResultsObj)))
            cat(sprintf("> variantType(res) <- c(%s)\n", paste(paste(names(variantType(vfResultsObj)), as.character(variantType(vfResultsObj)), sep="="), collapse=", ")))
        if (!all(variantLocation(vfResultsObj)))
            cat(sprintf("> variantLocation(res) <- c(%s)\n", paste(paste(names(variantLocation(vfResultsObj)), as.character(variantLocation(vfResultsObj)), sep="="), collapse=", ")))
        if (!all(variantConsequence(vfResultsObj)))
            cat(sprintf("> variantConsequence(res) <- c(%s)\n", paste(paste(names(variantConsequence(vfResultsObj)), as.character(variantConsequence(vfResultsObj)), sep="="), collapse=", ")))
        if (aaChangeType(vfResultsObj) != "Any")
          cat(sprintf("> aaChangeType(res) <- \"%s\"\n", aaChangeType(vfResultsObj)))
        if (minCUFC(vfResultsObj) > 0)
          cat(sprintf("> minCUFC(res) <- %.2f\n", minCUFC(vfResultsObj)))
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
        ## if (all(!is.na(crypsplice))) {
        ##   if (!is.na(minScore5ss(vfResultsObj)))
        ##     cat(sprintf("> minScore5ss(res) <- %.2f\n", minScore5ss(vfResultsObj)))
        ##   if (!is.na(minScore3ss(vfResultsObj)))
        ##     cat(sprintf("> minScore3ss(res) <- %.2f\n", minScore3ss(vfResultsObj)))
        ## }
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
