setMethod("show", signature(object="VariantFilteringResultsAIM"),
          function(object) {
            cat("\nVariantFiltering results object\n\n")
            cat(sprintf("  Variants segregate according to an %s inheritance model\n", inheritanceModel(object)))
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


setMethod("selectInheritancePattern", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@inheritancepattern
          })

setReplaceMethod("selectInheritancePattern", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value=c("None", "Autosomal Recessive Hom", "Autosomal Recessive Het",
                                     "Autosomal Dominant", "X-Linked", "de Novo")) {
                   value <- match.arg(value)
                   x@inheritancepattern <- as.character(value)
                   x
          })

setMethod("selectIndividual", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@indselected
          })

setReplaceMethod("selectIndividual", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                     x@indselected <- as.character(value)
                     x
          })

setMethod("selectIndexCase", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@selectindexcase
          })

setReplaceMethod("selectIndexCase", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                   x@selectindexcase <- as.character(value)
                   x
          })

                      ##############################
                  ### Autosomal Recessive Homozygous ###
                      ##############################

setMethod("selectCarrierRelative1", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@selectcarrierrelative1
          })

setReplaceMethod("selectCarrierRelative1", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                   x@selectcarrierrelative1 <- as.character(value)
                   x
          })

setMethod("selectCarrierRelative2", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@selectcarrierrelative2
          })

setReplaceMethod("selectCarrierRelative2", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                   x@selectcarrierrelative2 <- as.character(value)
                   x
          })

setMethod("selectAffectedRelative", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@selectaffectedrelative
          })

setReplaceMethod("selectAffectedRelative", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                   x@selectaffectedrelative <- as.character(value)
                   x
          })

                      ################################
                  ### Autosomal Recessive Heterozygous ###
                      ################################

setMethod("selectCarrierAllele1CH", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@selectcarrierallele1ch
          })

setReplaceMethod("selectCarrierAllele1CH", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                   x@selectcarrierallele1ch <- as.character(value)
                   x
          })

setMethod("selectCarrierAllele2CH", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@selectcarrierallele2ch
          })

setReplaceMethod("selectCarrierAllele2CH", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                   x@selectcarrierallele2ch <- as.character(value)
                   x
          })

setMethod("selectAffRelative1CH", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@selectaffrelative1ch
          })

setReplaceMethod("selectAffRelative1CH", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                   x@selectaffrelative1ch <- as.character(value)
                   x
          })


                      ##################
                  ### Autosomal Dominant ###
                      ##################

setMethod("selectUnaffectedRelative1AD", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@selectunaffectedrelative1ad
          })

setReplaceMethod("selectUnaffectedRelative1AD", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                   x@selectunaffectedrelative1ad <- as.character(value)
                   x
          })

setMethod("selectUnaffectedRelative2AD", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@selectunaffectedrelative2ad
          })

setReplaceMethod("selectUnaffectedRelative2AD", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                   x@selectunaffectedrelative2ad <- as.character(value)
                   x
          })

setMethod("selectAffectedRelative1AD", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@selectaffectedrelative1ad
          })

setReplaceMethod("selectAffectedRelative1AD", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                   x@selectaffectedrelative1ad <- as.character(value)
                   x
          })

                      ########
                  ### X-Linked ###
                      ########

setMethod("selectCarrRelFem1XL", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@selectcarrierrelativefemale1xl
          })

setReplaceMethod("selectCarrRelFem1XL", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                   x@selectcarrierrelativefemale1xl <- as.character(value)
                   x
          })

setMethod("selectAffRelMale1XL", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@selectaffectedrelativemale1xl
          })

setReplaceMethod("selectAffRelMale1XL", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                   x@selectaffectedrelativemale1xl <- as.character(value)
                   x
          })

setMethod("selectUnaffMale1XL", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@selectunaffectedmale1xl
          })

setReplaceMethod("selectUnaffMale1XL", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                   x@selectunaffectedmale1xl <- as.character(value)
                   x
          })

                      #######
                  ### de Novo ###
                      #######

setMethod("selectParent1DN", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@selectparent1dn
          })

setReplaceMethod("selectParent1DN", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                   x@selectparent1dn <- as.character(value)
                   x
          })

setMethod("selectParent2DN", signature(x="VariantFilteringResultsAIM"),
          function(x) {
              x@selectparent2dn
          })

setReplaceMethod("selectParent2DN", signature(x="VariantFilteringResultsAIM", value="ANY"),
                 function(x, value) {
                   x@selectparent2dn <- as.character(value)
                   x
          })


## get variants after applying all filters
setMethod("filteredVariants", signature(x="VariantFilteringResultsAIM"), 
          function(x, unusedColumns.rm=FALSE) {
            vars <- allVariants(x)
            rowsMask <- rep(TRUE, length(vars))

            if (selectInheritancePattern(x) == "Autosomal Recessive Hom") {
                maskIndexHomAlt <- sum(which(vars$HomozygousAlt == selectIndexCase(x))) > 0
                autosomalChrVars <- seqnames(vars) != "chrX" & seqnames(vars) != "chrY"
                rowsMask <- rowsMask & maskIndexHomAlt & autosomalChrVars
                if (!is.na(selectCarrierRelative1(x))) {
                    maskCarrier1Het <- sum(which(vars$Heterozygous == selectCarrierRelative1(x))) > 0
                    rowsMask <- rowsMask & maskCarrier1Het
                }
                if (!is.na(selectCarrierRelative2(x))) {
                    maskCarrier2Het <- sum(which(vars$Heterozygous == selectCarrierRelative2(x))) > 0
                    rowsMask <- rowsMask & maskCarrier2Het
                }
                if (!is.na(selectAffectedRelative(x))) {
                    maskRelativeHomAlt <- sum(which(vars$HomozygousAlt == selectAffectedRelative(x))) > 0
                    rowsMask <- rowsMask & maskRelativeHomAlt
                }
            }

            if (selectInheritancePattern(x) == "Autosomal Recessive Het") {
                maskIndexHet <- sum(which(vars$Heterozygous == selectIndexCase(x))) > 0
                maskCoding <- vars$LOCATION == "coding"
                maskIndexHetCod <- maskIndexHet & maskCoding
                maskAllele1Het <- sum(which(vars$Heterozygous == selectCarrierAllele1CH(x))) > 0
                maskAllele2Het <- sum(which(vars$Heterozygous == selectCarrierAllele2CH(x))) > 0

                ## maskIndexHet <- sum(which(vars$Heterozygous == "NA12878")) > 0
                ## maskCoding <- vars$LOCATION == "coding"
                ## maskIndexHetCod <- maskIndexHet & maskCoding
                ## maskAllele1Het <- sum(which(vars$Heterozygous == "NA12891")) > 0
                ## maskAllele2Het <- sum(which(vars$Heterozygous == "NA12892")) > 0

                ## pull out common variants between the two carriers
                maskCommonHet <- maskAllele1Het & maskAllele2Het
                maskNotCommon1 <- maskAllele1Het & !maskCommonHet
                maskNotCommon2 <- maskAllele2Het & !maskCommonHet

                maskIndexWith1 <- maskIndexHetCod & maskNotCommon1
                maskIndexWith2 <- maskIndexHetCod & maskNotCommon2
                
                ## keep genes which appears in both individuals once pulled out common variants (i.e. same genes affected by different variants)
                genes <- unique(vars$GENE[maskIndexWith1])[unique(vars$GENE[maskIndexWith1]) %in% unique(vars$GENE[maskIndexWith2])]
                maskCommonGenes <- vars$GENE %in% genes
                
                maskIndex12 <- maskIndexWith1 | maskIndexWith2
                maskCompHet <- maskIndex12 & maskCommonGenes

                rowsMask <- rowsMask & maskCompHet
                if (!is.na(selectAffRelative1CH(x))) {
                    maskAffRel1Het <- sum(which(vars$Heterozygous == selectAffRelative1CH(x))) > 0
                    rowsMask <- rowsMask & maskAffRel1Het
                }
            }
            

            if (selectInheritancePattern(x) == "Autosomal Dominant") {
                maskIndexHet <- sum(which(vars$Heterozygous == selectIndexCase(x))) > 0
                autosomalChrVars <- seqnames(vars) != "chrX" & seqnames(vars) != "chrY"                
                rowsMask <- rowsMask & maskIndexHet & autosomalChrVars
                if (!is.na(selectUnaffectedRelative1AD(x))) {
                    maskUnaff1HomRef <- sum(which(vars$HomozygousRef == selectUnaffectedRelative1AD(x))) > 0
                    rowsMask <- rowsMask & maskUnaff1HomRef
                }
                if (!is.na(selectUnaffectedRelative2AD(x))) {
                    maskUnaff2HomRef <- sum(which(vars$HomozygousRef == selectUnaffectedRelative2AD(x))) > 0
                    rowsMask <- rowsMask & maskUnaff2HomRef
                }
                if (!is.na(selectAffectedRelative1AD(x))) {
                    maskRelative1Het <- sum(which(vars$Heterozygous == selectAffectedRelative1AD(x))) > 0
                    rowsMask <- rowsMask & maskRelative1Het
                }
            }

            if (selectInheritancePattern(x) == "X-Linked") {
                maskIndexHomAlt <- sum(which(vars$HomozygousAlt == selectIndexCase(x))) > 0
                xChrVars <- seqnames(vars) == "chrX"
                rowsMask <- rowsMask & maskIndexHomAlt & xChrVars
                if (!is.na(selectCarrRelFem1XL(x))) {
                    maskCarrierFemale1Het <- sum(which(vars$Heterozygous == selectCarrRelFem1XL(x))) > 0
                    rowsMask <- rowsMask & maskCarrierFemale1Het
                }
                if (!is.na(selectAffRelMale1XL(x))) {
                    maskAffRelMale1HomAlt <- sum(which(vars$HomozygousAlt == selectAffRelMale1XL(x))) > 0
                    rowsMask <- rowsMask & maskAffRelMale1HomAlt
                }
                if (!is.na(selectUnaffMale1XL(x))) {
                    maskUnaffMale1HomRef <- sum(which(vars$HomozygousRef == selectUnaffMale1XL(x))) > 0
                    rowsMask <- rowsMask & maskUnaffMale1HomRef
                }
            }

            if (selectInheritancePattern(x) == "de Novo") {
                ## here we can not just do != HomRef because that would also include the "no calls"
                ## we expect that a de novo variant will be shown as heterozygous in the autosomals
                ## and x chromosomes in females and homozygous alt in the x chr in males
                ## the user has to be aware about the sex of the index individual 
                maskIndexHet <- sum(which(vars$Heterozygous == selectIndexCase(x))) > 0
                maskIndexHomAlt <- sum(which(vars$HomozygousAlt == selectIndexCase(x))) > 0                
                xChrVars <- seqnames(vars) == "chrX"
                maskXChrHomAlt <- maskIndexHomAlt & xChrVars
                maskDeNovoVars <- maskIndexHet | maskXChrHomAlt
                rowsMask <- rowsMask & maskDeNovoVars
                if (!is.na(selectParent1DN(x))) {
                    maskParent1HomRef <- sum(which(vars$HomozygousRef == selectParent1DN(x))) > 0
                    rowsMask <- rowsMask & maskParent1HomRef
                }
                if (!is.na(selectParent2DN(x))) {
                    maskParent2HomRef <- sum(which(vars$HomozygousRef == selectParent2DN(x))) > 0
                    rowsMask <- rowsMask & maskParent2HomRef
                }
            }                
            

            ## individual
            ## 
            if (!is.na(selectIndividual(x))) {
                maskIndHom <- sum(which(vars$HomozygousAlt == selectIndividual(x))) > 0
                maskIndHet <- sum(which(vars$Heterozygous == selectIndividual(x))) > 0
                rowsMask <- rowsMask & (maskIndHet | maskIndHom)
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
 

            ## location
            
            ##
            ## filter by location
            ##

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
                if (!namaf(X))
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


setMethod("reportVariants", signature(vfResultsObj="VariantFilteringResultsAIM"),
          function(vfResultsObj, type=c("shiny", "csv", "tsv"), file=NULL) {
              
  type <- match.arg(type)

  if (class(vfResultsObj) != "VariantFilteringResultsAIM")
    stop("Input argument 'vfResultsObj' should be an object of class 'VariantFilteringResultsAIM'.")

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

                ## inheritance tab
                conditionalPanel(condition="input.tsp == 'inheritance'",
                                 selectInput("selectInheritancePattern", "Inheritance Model:",
                                             choices=c("None", "Autosomal Recessive Hom", "Autosomal Recessive Het",
                                     "Autosomal Dominant", "X-Linked", "de Novo"))),

                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern != 'None'",
                                 selectInput("selectIndexCase", "Index case:",
                                             choices=indlist)),
                  
                      ##############################
                  ### Autosomal Recessive Homozygous ###
                      ##############################
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Recessive Hom'",
                                 checkboxInput("carrierRelative1Flag",
                                               "Carrier relative 1", FALSE)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Recessive Hom' &&
                                            input.carrierRelative1Flag == true",
                                 selectInput("selectCarrierRelative1", "Carrier relative 1:",
                                             choices=indlist)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Recessive Hom'",
                                 checkboxInput("carrierRelative2Flag",
                                               "Carrier relative 2", FALSE)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Recessive Hom' &&
                                            input.carrierRelative2Flag == true",
                                 selectInput("selectCarrierRelative2", "Carrier relative 2:",
                                             choices=indlist)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Recessive Hom'",
                                 checkboxInput("affectedRelativeFlag",
                                               "Affected relative", FALSE)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Recessive Hom' &&
                                            input.affectedRelativeFlag == true",
                                 selectInput("selectAffectedRelative", "Affected relative:",
                                             choices=indlist)),
                                               
                      ################################
                  ### Autosomal Recessive Heterozygous ###
                      ################################

                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Recessive Het'",
                                 checkboxInput("carrierAllele1FlagCH",
                                               "Carrier allele 1", FALSE)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Recessive Het' &&
                                            input.carrierAllele1FlagCH == true",
                                 selectInput("selectCarrierAllele1CH", "Carrier allele 1:",
                                             choices=indlist)),

                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Recessive Het'",
                                 checkboxInput("carrierAllele2FlagCH",
                                               "Carrier allele 2", FALSE)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Recessive Het' &&
                                            input.carrierAllele2FlagCH == true",
                                 selectInput("selectCarrierAllele2CH", "Carrier allele 2:",
                                             choices=indlist)),

                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Recessive Het'",
                                 checkboxInput("affectedRelative1FlagCH",
                                               "Affected relative 1", FALSE)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Recessive Het' &&
                                            input.affectedRelative1FlagCH == true",
                                 selectInput("selectAffRelative1CH", "Affected relative 1:",
                                             choices=indlist)),



                      ##################
                  ### Autosomal Dominant ###
                      ##################

                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Dominant'",
                                 checkboxInput("unaffectedRelative1FlagAD",
                                               "Unaffected relative 1", FALSE)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Dominant' &&
                                            input.unaffectedRelative1FlagAD == true",
                                 selectInput("selectUnaffectedRelative1AD", "Unaffected relative 1:",
                                             choices=indlist)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Dominant'",
                                 checkboxInput("unaffectedRelative2FlagAD",
                                               "Unaffected relative 2", FALSE)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Dominant' &&
                                            input.unaffectedRelative2FlagAD == true",
                                 selectInput("selectUnaffectedRelative2AD", "Unaffected relative 2:",
                                             choices=indlist)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Dominant'",
                                 checkboxInput("affectedRelativeFlag1AD",
                                               "Affected relative 1", FALSE)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'Autosomal Dominant' &&
                                            input.affectedRelativeFlag1AD == true",
                                 selectInput("selectAffectedRelative1AD", "Affected relative:",
                                             choices=indlist)),
                  
                      ########
                  ### X-Linked ###
                      ########

                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'X-Linked'",
                                 checkboxInput("carrierRelativeFemale1FlagXL",
                                               "Carrier relative female 1", FALSE)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'X-Linked' &&
                                            input.carrierRelativeFemale1FlagXL == true",
                                 selectInput("selectCarrRelFem1XL", "Carrier relative female 1:",
                                             choices=indlist)),

                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'X-Linked'",
                                 checkboxInput("affectedRelativeMale1FlagXL",
                                               "Affected relative male 1", FALSE)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'X-Linked' &&
                                            input.affectedRelativeMale1FlagXL == true",
                                 selectInput("selectAffRelMale1XL", "Affected relative male 1:",
                                             choices=indlist)),
                
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'X-Linked'",
                                 checkboxInput("unaffectedMale1FlagXL",
                                               "Unaffected male 1", FALSE)),
                  
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'X-Linked' &&
                                            input.unaffectedMale1FlagXL == true",
                                 selectInput("selectUnaffMale1XL", "Unaffected male 1:",
                                             choices=indlist)),
                  
                      #######
                  ### de Novo ###
                      #######

                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'de Novo'",
                                 checkboxInput("parent1FlagDN",
                                               "Parent 1", FALSE)),
                   
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'de Novo' &&
                                            input.parent1FlagDN == true",
                                 selectInput("selectParent1DN", "Parent 1:",
                                             choices=indlist)),
                   
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'de Novo'",
                                 checkboxInput("parent2FlagDN",
                                               "Parent 2", FALSE)),
                   
                conditionalPanel(condition="input.tsp == 'inheritance' &&
                                            input.selectInheritancePattern == 'de Novo' &&
                                            input.parent2FlagDN == true",
                                 selectInput("selectParent2DN", "Parent 2:",
                                             choices=indlist)),
                    
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
                  tabPanel("Inheritance", tableOutput('tableInheritance'), value="inheritance"), 
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
          
        ## inheritance pattern
        selectInheritancePattern(vfResultsObj) <- input$selectInheritancePattern

        ## index case
        selectIndexCase(vfResultsObj) <- input$selectIndexCase

                      ##############################
                  ### Autosomal Recessive Homozygous ###
                      ##############################
        
        ## carrier relative 1
          if (input$carrierRelative1Flag)
              selectCarrierRelative1(vfResultsObj) <- input$selectCarrierRelative1
          else
              selectCarrierRelative1(vfResultsObj) <- NA_character_
        
        ## carrier relative 2
          if (input$carrierRelative2Flag)
              selectCarrierRelative2(vfResultsObj) <- input$selectCarrierRelative2
          else
              selectCarrierRelative2(vfResultsObj) <- NA_character_
        
        ## affected relative
          if (input$affectedRelativeFlag)
              selectAffectedRelative(vfResultsObj) <- input$selectAffectedRelative
          else
              selectAffectedRelative(vfResultsObj) <- NA_character_

        
                      ################################
                  ### Autosomal Recessive Heterozygous ###
                      ################################

        ## carrier allele 1
          if (input$carrierAllele1FlagCH)
              selectCarrierAllele1CH(vfResultsObj) <- input$selectCarrierAllele1CH
          else
              selectCarrierAllele1CH(vfResultsObj) <- NA_character_
 
        ## carrier allele 2
          if (input$carrierAllele2FlagCH)
              selectCarrierAllele2CH(vfResultsObj) <- input$selectCarrierAllele2CH
          else
              selectCarrierAllele2CH(vfResultsObj) <- NA_character_
 
        ## affected relative 1
          if (input$affectedRelative1FlagCH)
              selectAffRelative1CH(vfResultsObj) <- input$selectAffRelative1CH
          else
              selectAffRelative1CH(vfResultsObj) <- NA_character_

        
                      ##################
                  ### Autosomal Dominant ###
                      ##################
        
        ## unaffected relative 1
          if (input$unaffectedRelative1FlagAD)
              selectUnaffectedRelative1AD(vfResultsObj) <- input$selectUnaffectedRelative1AD
          else
              selectUnaffectedRelative1AD(vfResultsObj) <- NA_character_
        
        ## unaffected relative 2
          if (input$unaffectedRelative2FlagAD)
              selectUnaffectedRelative2AD(vfResultsObj) <- input$selectUnaffectedRelative2AD
          else
              selectUnaffectedRelative2AD(vfResultsObj) <- NA_character_
        
        ## affected relative 1
          if (input$affectedRelativeFlag1AD)
              selectAffectedRelative1AD(vfResultsObj) <- input$selectAffectedRelative1AD
          else
              selectAffectedRelative1AD(vfResultsObj) <- NA_character_

        
                      ########
                  ### X-Linked ###
                      ########

        ## carrier relative female 1
          if (input$carrierRelativeFemale1FlagXL)
              selectCarrRelFem1XL(vfResultsObj) <- input$selectCarrRelFem1XL
          else
              selectCarrRelFem1XL(vfResultsObj) <- NA_character_

        ## affected relative male 1
          if (input$affectedRelativeMale1FlagXL)
              selectAffRelMale1XL(vfResultsObj) <- input$selectAffRelMale1XL
          else
              selectAffRelMale1XL(vfResultsObj) <- NA_character_

        ## unaffected male 1
          if (input$unaffectedMale1FlagXL)
              selectUnaffMale1XL(vfResultsObj) <- input$selectUnaffMale1XL
          else
              selectUnaffMale1XL(vfResultsObj) <- NA_character_

        
                      #######
                  ### de Novo ###
                      #######
        
        ## parent 1
          if (input$parent1FlagDN)
              selectParent1DN(vfResultsObj) <- input$selectParent1DN
          else
              selectParent1DN(vfResultsObj) <- NA_character_

        ## parent 2
          if (input$parent2FlagDN)
              selectParent2DN(vfResultsObj) <- input$selectParent2DN
          else
              selectParent2DN(vfResultsObj) <- NA_character_
        
      ########################################
        
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
        
        ## inheritance pattern
        selectInheritancePattern(vfResultsObj) <- input$selectInheritancePattern
        
        ## index case
        selectIndexCase(vfResultsObj) <- input$selectIndexCase

                      ##############################
                  ### Autosomal Recessive Homozygous ###
                      ##############################
        
        ## carrier relative 1
          if (input$carrierRelative1Flag)
              selectCarrierRelative1(vfResultsObj) <- input$selectCarrierRelative1
          else
              selectCarrierRelative1(vfResultsObj) <- NA_character_
        
        ## carrier relative 2
          if (input$carrierRelative2Flag)
              selectCarrierRelative2(vfResultsObj) <- input$selectCarrierRelative2
          else
              selectCarrierRelative2(vfResultsObj) <- NA_character_
        
        ## affected relative
          if (input$affectedRelativeFlag)
              selectAffectedRelative(vfResultsObj) <- input$selectAffectedRelative
          else
              selectAffectedRelative(vfResultsObj) <- NA_character_

                      ################################
                  ### Autosomal Recessive Heterozygous ###
                      ################################

        ## carrier allele 1
          if (input$carrierAllele1FlagCH)
              selectCarrierAllele1CH(vfResultsObj) <- input$selectCarrierAllele1CH
          else
              selectCarrierAllele1CH(vfResultsObj) <- NA_character_
 
        ## carrier allele 2
          if (input$carrierAllele2FlagCH)
              selectCarrierAllele2CH(vfResultsObj) <- input$selectCarrierAllele2CH
          else
              selectCarrierAllele2CH(vfResultsObj) <- NA_character_
 
        ## affected relative 1
          if (input$affectedRelative1FlagCH)
              selectAffRelative1CH(vfResultsObj) <- input$selectAffRelative1CH
          else
              selectAffRelative1CH(vfResultsObj) <- NA_character_


                      ##################
                  ### Autosomal Dominant ###
                      ##################
        
        ## unaffected relative 1
          if (input$unaffectedRelative1FlagAD)
              selectUnaffectedRelative1AD(vfResultsObj) <- input$selectUnaffectedRelative1AD
          else
              selectUnaffectedRelative1AD(vfResultsObj) <- NA_character_
        
        ## unaffected relative 2
          if (input$unaffectedRelative2FlagAD)
              selectUnaffectedRelative2AD(vfResultsObj) <- input$selectUnaffectedRelative2AD
          else
              selectUnaffectedRelative2AD(vfResultsObj) <- NA_character_
        
        ## affected relative 1
          if (input$affectedRelativeFlag1AD)
              selectAffectedRelative1AD(vfResultsObj) <- input$selectAffectedRelative1AD
          else
              selectAffectedRelative1AD(vfResultsObj) <- NA_character_

                      ########
                  ### X-Linked ###
                      ########

        ## carrier relative female 1
          if (input$carrierRelativeFemale1FlagXL)
              selectCarrRelFem1XL(vfResultsObj) <- input$selectCarrRelFem1XL
          else
              selectCarrRelFem1XL(vfResultsObj) <- NA_character_

        ## affected relative male 1
          if (input$affectedRelativeMale1FlagXL)
              selectAffRelMale1XL(vfResultsObj) <- input$selectAffRelMale1XL
          else
              selectAffRelMale1XL(vfResultsObj) <- NA_character_

        ## unaffected male 1
          if (input$unaffectedMale1FlagXL)
              selectUnaffMale1XL(vfResultsObj) <- input$selectUnaffMale1XL
          else
              selectUnaffMale1XL(vfResultsObj) <- NA_character_
        
                      #######
                  ### de Novo ###
                      #######
        
        ## parent 1
          if (input$parent1FlagDN)
              selectParent1DN(vfResultsObj) <- input$selectParent1DN
          else
              selectParent1DN(vfResultsObj) <- NA_character_

        ## parent 2
          if (input$parent2FlagDN)
              selectParent2DN(vfResultsObj) <- input$selectParent2DN
          else
              selectParent2DN(vfResultsObj) <- NA_character_

      ########################################
         
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
      vdf$HomozygousRef <- sapply(vdf$HomozygousRef, paste, collapse=", ")
     
      
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

    output$tableInheritance <- renderTable({
        filteredVariantsReact()[, c("POSITION", "CDS", "HomozygousRef", "Heterozygous", "HomozygousAlt")]
    }, NA.string="NA", sanitize.text.function=function(x){x})

    output$tableGenome <- renderTable({
      filteredVariantsReact()[, c("VarID", "POSITION", "dbSNP", "TYPE")]
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
          
        ## inheritance pattern
        selectInheritancePattern(vfResultsObj) <- input$selectInheritancePattern

        ## index case
        selectIndexCase(vfResultsObj) <- input$selectIndexCase
        
                      ##############################
                  ### Autosomal Recessive Homozygous ###
                      ##############################
        
        ## carrier relative 1
          if (input$carrierRelative1Flag)
              selectCarrierRelative1(vfResultsObj) <- input$selectCarrierRelative1
          else
              selectCarrierRelative1(vfResultsObj) <- NA_character_
        ## carrier relative 2
          if (input$carrierRelative2Flag)
              selectCarrierRelative2(vfResultsObj) <- input$selectCarrierRelative2
          else
              selectCarrierRelative2(vfResultsObj) <- NA_character_
        ## affected relative
          if (input$affectedRelativeFlag)
              selectAffectedRelative(vfResultsObj) <- input$selectAffectedRelative
          else
              selectAffectedRelative(vfResultsObj) <- NA_character_

                      ################################
                  ### Autosomal Recessive Heterozygous ###
                      ################################

        ## carrier allele 1
          if (input$carrierAllele1FlagCH)
              selectCarrierAllele1CH(vfResultsObj) <- input$selectCarrierAllele1CH
          else
              selectCarrierAllele1CH(vfResultsObj) <- NA_character_
 
        ## carrier allele 2
          if (input$carrierAllele2FlagCH)
              selectCarrierAllele2CH(vfResultsObj) <- input$selectCarrierAllele2CH
          else
              selectCarrierAllele2CH(vfResultsObj) <- NA_character_
 
        ## affected relative 1
          if (input$affectedRelative1FlagCH)
              selectAffRelative1CH(vfResultsObj) <- input$selectAffRelative1CH
          else
              selectAffRelative1CH(vfResultsObj) <- NA_character_


                      ##################
                  ### Autosomal Dominant ###
                      ##################
        
        ## unaffected relative 1
          if (input$unaffectedRelative1FlagAD)
              selectUnaffectedRelative1AD(vfResultsObj) <- input$selectUnaffectedRelative1AD
          else
              selectUnaffectedRelative1AD(vfResultsObj) <- NA_character_
        
        ## unaffected relative 2
          if (input$unaffectedRelative2FlagAD)
              selectUnaffectedRelative2AD(vfResultsObj) <- input$selectUnaffectedRelative2AD
          else
              selectUnaffectedRelative2AD(vfResultsObj) <- NA_character_
        
        ## affected relative 1
          if (input$affectedRelativeFlag1AD)
              selectAffectedRelative1AD(vfResultsObj) <- input$selectAffectedRelative1AD
          else
              selectAffectedRelative1AD(vfResultsObj) <- NA_character_

                      ########
                  ### X-Linked ###
                      ########

        ## carrier relative female 1
          if (input$carrierRelativeFemale1FlagXL)
              selectCarrRelFem1XL(vfResultsObj) <- input$selectCarrRelFem1XL
          else
              selectCarrRelFem1XL(vfResultsObj) <- NA_character_

        ## affected relative male 1
          if (input$affectedRelativeMale1FlagXL)
              selectAffRelMale1XL(vfResultsObj) <- input$selectAffRelMale1XL
          else
              selectAffRelMale1XL(vfResultsObj) <- NA_character_

        ## unaffected male 1
          if (input$unaffectedMale1FlagXL)
              selectUnaffMale1XL(vfResultsObj) <- input$selectUnaffMale1XL
          else
              selectUnaffMale1XL(vfResultsObj) <- NA_character_
        
                      #######
                  ### de Novo ###
                      #######
        
        ## parent 1
          if (input$parent1FlagDN)
              selectParent1DN(vfResultsObj) <- input$selectParent1DN
          else
              selectParent1DN(vfResultsObj) <- NA_character_

        ## parent 2
          if (input$parent2FlagDN)
              selectParent2DN(vfResultsObj) <- input$selectParent2DN
          else
              selectParent2DN(vfResultsObj) <- NA_character_

      ########################################
        
 
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
          
        ## inheritance pattern
        selectInheritancePattern(vfResultsObj) <- input$selectInheritancePattern
        
        ## index case
        selectIndexCase(vfResultsObj) <- input$selectIndexCase
        
                      ##############################
                  ### Autosomal Recessive Homozygous ###
                      ##############################

        ## carrier relative 1
          if (input$carrierRelative1Flag)
              selectCarrierRelative1(vfResultsObj) <- input$selectCarrierRelative1
          else
              selectCarrierRelative1(vfResultsObj) <- NA_character_
        
        ## carrier relative 2
          if (input$carrierRelative2Flag)
              selectCarrierRelative2(vfResultsObj) <- input$selectCarrierRelative2
          else
              selectCarrierRelative2(vfResultsObj) <- NA_character_
        
        ## affected relative
          if (input$affectedRelativeFlag)
              selectAffectedRelative(vfResultsObj) <- input$selectAffectedRelative
          else
              selectAffectedRelative(vfResultsObj) <- NA_character_
        
                      ################################
                  ### Autosomal Recessive Heterozygous ###
                      ################################

        ## carrier allele 1
          if (input$carrierAllele1FlagCH)
              selectCarrierAllele1CH(vfResultsObj) <- input$selectCarrierAllele1CH
          else
              selectCarrierAllele1CH(vfResultsObj) <- NA_character_
 
        ## carrier allele 2
          if (input$carrierAllele2FlagCH)
              selectCarrierAllele2CH(vfResultsObj) <- input$selectCarrierAllele2CH
          else
              selectCarrierAllele2CH(vfResultsObj) <- NA_character_
 
        ## affected relative 1
          if (input$affectedRelative1FlagCH)
              selectAffRelative1CH(vfResultsObj) <- input$selectAffRelative1CH
          else
              selectAffRelative1CH(vfResultsObj) <- NA_character_


                      ##################
                  ### Autosomal Dominant ###
                      ##################
        
        ## unaffected relative 1
          if (input$unaffectedRelative1FlagAD)
              selectUnaffectedRelative1AD(vfResultsObj) <- input$selectUnaffectedRelative1AD
          else
              selectUnaffectedRelative1AD(vfResultsObj) <- NA_character_
        
        ## unaffected relative 2
          if (input$unaffectedRelative2FlagAD)
              selectUnaffectedRelative2AD(vfResultsObj) <- input$selectUnaffectedRelative2AD
          else
              selectUnaffectedRelative2AD(vfResultsObj) <- NA_character_
        
        ## affected relative 1
          if (input$affectedRelativeFlag1AD)
              selectAffectedRelative1AD(vfResultsObj) <- input$selectAffectedRelative1AD
          else
              selectAffectedRelative1AD(vfResultsObj) <- NA_character_

                      ########
                  ### X-Linked ###
                      ########

        ## carrier relative female 1
          if (input$carrierRelativeFemale1FlagXL)
              selectCarrRelFem1XL(vfResultsObj) <- input$selectCarrRelFem1XL
          else
              selectCarrRelFem1XL(vfResultsObj) <- NA_character_

        ## affected relative male 1
          if (input$affectedRelativeMale1FlagXL)
              selectAffRelMale1XL(vfResultsObj) <- input$selectAffRelMale1XL
          else
              selectAffRelMale1XL(vfResultsObj) <- NA_character_

        ## unaffected male 1
          if (input$unaffectedMale1FlagXL)
              selectUnaffMale1XL(vfResultsObj) <- input$selectUnaffMale1XL
          else
              selectUnaffMale1XL(vfResultsObj) <- NA_character_
        
                      #######
                  ### de Novo ###
                      #######
        
        ## parent 1
          if (input$parent1FlagDN)
              selectParent1DN(vfResultsObj) <- input$selectParent1DN
          else
              selectParent1DN(vfResultsObj) <- NA_character_

        ## parent 2
          if (input$parent2FlagDN)
              selectParent2DN(vfResultsObj) <- input$selectParent2DN
          else
              selectParent2DN(vfResultsObj) <- NA_character_

      ########################################
 
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
