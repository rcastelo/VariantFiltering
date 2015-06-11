## main function carrying out the annotation of genetic variants taken as an
## input GRanges object. It uses different functions from this package and
## from the VariantAnnotation package

annotationEngine <- function(variantsVR, param, cache=new.env(parent=emptyenv()),
                             BPPARAM=bpparam("SerialParam")) {
  
  if (length(variantsVR) == 0) {
    variantsVR_annotated <- variantsVR
    mcols(variantsVR_annotated) <- .emptyAnnotations()
    return(variantsVR_annotated)
  }

  bsgenome <- param$bsgenome
  orgdb <- param$orgdb
  txdb <- param$txdb
  snpdb <- param$snpdb
  spliceSiteMatrices <- param$spliceSiteMatrices
  radicalAAchangeMatrix <- param$radicalAAchangeMatrix
  allTranscripts <- param$allTranscripts
  otherAnnotations <- param$otherAnnotations
  codonusageTable <- param$codonusageTable
  geneticCode <- param$geneticCode

  ##############################
  ##                          ##
  ## CLEAN UP VARIANT INFO    ##
  ##                          ##
  ##############################

  ## clean the variant information landed on the 'VARID' slot of the input variantsVR 'GenomicRanges'
  ## object to leave either dbSNP ids given to the input VCF in this field or a maximum of 20
  ## characters. When the dbSNP annotation does not succeed in fetching information, this 'VARID'
  ## coming from the 'ID' column of the original VCF file will be used for that purpose

  vnames <- variantsVR$VARID
  if (!is.null(vnames)) {
    vnames2 <- vnames
    vnames2 <- rep(NA_character_, length(vnames))

    mt <- gregexpr("^[a-zA-Z0-9]+:[0-9]+_[ACGT/]+", vnames)
    mtstart <- unlist(mt, use.names=FALSE)
    mtlength <- sapply(mt, attr, "match.length")
    vnames2[mtstart != -1] <- substr(vnames[mtstart != -1], mtstart[mtstart != -1], mtlength[mtstart != -1])

    mt <- gregexpr("RS=[a-zA-Z]+[0-9]+", vnames)
    mtstart <- sapply(mt, "[", 1)
    mtlength <- sapply(mt, attr, "match.length")
    mtlength <- sapply(mtlength, "[", 1)
    vnames2[mtstart != -1] <- substr(vnames[mtstart != -1], mtstart[mtstart != -1]+3, mtlength[mtstart != -1])

    mt <- gregexpr("[a-zA-Z]+[0-9]+", vnames)
    mtstart <- sapply(mt, "[", 1)
    mtlength <- sapply(mt, attr, "match.length")
    mtlength <- sapply(mtlength, "[", 1)
    vnames2[mtstart != -1] <- substr(vnames[mtstart != -1], mtstart[mtstart != -1], mtlength[mtstart != -1])

    wh <- nchar(vnames2) > 20
    vnames2[wh] <- paste0(substr(vnames2[wh], 1, 20), "...")

    variantsVR$VARID <- vnames2
  }

  ##############################
  ##                          ##
  ## ANNOTATE TYPE OF VARIANT ##
  ##                          ##
  ##############################
  
  message("Annotating variant type (SNV, Insertion, Deletion, MNV, Delins)")
  mcols(variantsVR) <- cbind(mcols(variantsVR), typeOfVariants(variantsVR))

  ##############################
  ##                          ##
  ##  SNP-CENTRIC ANNOTATIONS ##
  ##                          ##
  ##############################

  ## do this before variants get replicated because of different functional annotations
  variantsVR$dbSNP <- rep(NA_character_, times=length(variantsVR))
  for (i in seq_len(length(snpdb))) {
    message(sprintf("Annotating dbSNP identifiers with %s", names(snpdb)[i]))
    res <- annotateVariants(snpdb[[i]], variantsVR, param, BPPARAM=BPPARAM)
    maskNAdbsnp <- is.na(res$dbSNP)
    maskNAannotdbsnp <- is.na(variantsVR$dbSNP)
    variantsVR$dbSNP[maskNAannotdbsnp & !maskNAdbsnp] <- res$dbSNP[maskNAannotdbsnp & !maskNAdbsnp]
    variantsVR$dbSNP[!maskNAannotdbsnp & !maskNAdbsnp] <-
      paste(variantsVR$dbSNP[!maskNAannotdbsnp & !maskNAdbsnp], res$dbSNP[!maskNAannotdbsnp & !maskNAdbsnp], sep=", ")
  }

  #######################
  ##                   ##
  ## ANNOTATE LOCATION ##
  ##                   ##
  #######################

  ## at the moment we are not interested in intergenic variants and we also leave promoter region
  ## boundaries at their default value. This could be parametrized if needed by the 'VariantFilteringParam' input object
  message("Annotating location with VariantAnnotation::locateVariants()")
  located_variantsVR <- .locateAllVariants(vfParam=param, query=as(variantsVR, "GRanges"),
                                           subject=txdb, cache=cache, BPPARAM=BPPARAM)
  ## located_variantsVR <- locateVariants(query=as(variantsVR, "GRanges"), subject=txdb,
  ##                                      region=AllVariants(intergenic=IntergenicVariants(0, 0)),
  ##                                      cache=cache)
  variantsVR_annotated <- variantsVR[located_variantsVR$QUERYID] ## REPLACE variantsVR_annotated by variantsVR ???
  variantsVR_annotated$LOCATION <- located_variantsVR$LOCATION
  variantsVR_annotated$LOCSTART <- located_variantsVR$LOCSTART
  variantsVR_annotated$LOCEND <- located_variantsVR$LOCEND
  variantsVR_annotated$LOCSTRAND <- strand(located_variantsVR)
  variantsVR_annotated$QUERYID <- located_variantsVR$QUERYID
  variantsVR_annotated$TXID <- located_variantsVR$TXID
  variantsVR_annotated$CDSID <- located_variantsVR$CDSID
  variantsVR_annotated$GENEID <- located_variantsVR$GENEID
  rm(located_variantsVR)

  ## if the argument 'allTranscripts' is set to 'FALSE' then keep only once identical variants
  ## annotated with the same type of region from the same gene but in different tx
  if (!allTranscripts) {
    selmcols <- c("QUERYID", "LOCATION", "GENEID")
    dupsmask <- duplicated(mcols(variantsVR_annotated)[, selmcols])
    variantsVR_annotated <- variantsVR_annotated[!dupsmask]
  }

  rmcols <- match("QUERYID", colnames(mcols(variantsVR_annotated)))
  mcols(variantsVR_annotated) <- mcols(variantsVR_annotated)[, -rmcols]

  ## annotate cDNA position where applicable
  maskexonic <- variantsVR_annotated$LOCATION %in% c("coding", "fiveUTR", "threeUTR")
  cDNAloc <- .cDNAloc(variantsVR_annotated[maskexonic], txdb)
  variantsVR_annotated$cDNALOC <- rep(NA_integer_, length(variantsVR_annotated))
  variantsVR_annotated$cDNALOC[maskexonic] <- cDNAloc

  ##############################
  ##                          ##
  ## ANNOTATE CODING VARIANTS ##
  ##                          ##
  ##############################
  
  message("Annotating coding variants VariantAnnotation::predictCoding()")
  ## following the help page of 'predictCoding()', expand to have one row per alternate allele
  variantsVR_annotated_coding <- VRanges()
  if (any(variantsVR_annotated$LOCATION == "coding")) {
    variantsVR_annotated_coding <- variantsVR_annotated[variantsVR_annotated$LOCATION == "coding"]

    ## set temporarily the sequence style to the one of the genome package to
    ## access the genome sequence without chromosome nomenclature problems
    origVarLevelsStyle <- seqlevelsStyle(variantsVR_annotated_coding)
    seqlevelsStyle(variantsVR_annotated_coding) <- seqlevelsStyle(bsgenome)
    origTxDbLevelsStyle <- seqlevelsStyle(txdb)
    seqlevelsStyle(txdb) <- seqlevelsStyle(bsgenome)

    commonChr <- intersect(seqlevels(variantsVR_annotated_coding), seqlevels(bsgenome))
    if (any(genome(variantsVR_annotated_coding)[commonChr] != genome(bsgenome)[commonChr])) {
      warning(sprintf("Assumming %s represent the same genome build.",
                      paste(c(unique(genome(variantsVR_annotated_coding)[commonChr]), unique(genome(bsgenome)[commonChr])),
                            collapse=" and ")))
      ## this generates conflicts with genome(txdb) in predictCoding() but it works
      ## we do not do anything about it and we just provide the warning
      ## genome(variantsVR_annotated_coding_exp) <- genome(bsgenome)
    }

    variantsVR_annotated_coding_exp <- as(variantsVR_annotated_coding, "GRanges")
    rmcols <- match(c("TXID", "CDSID", "GENEID"), colnames(mcols(variantsVR_annotated_coding_exp)))
    mcols(variantsVR_annotated_coding_exp) <- mcols(variantsVR_annotated_coding_exp)[, -rmcols]

    GRanges_coding_uq <- predictCoding(query=variantsVR_annotated_coding_exp,
                                       subject=txdb, seqSource=bsgenome, genetic.code=geneticCode,
                                       varAllele=DNAStringSet(alt(variantsVR_annotated_coding)))

    ## set back the sequence style
    seqlevelsStyle(variantsVR_annotated_coding) <- origVarLevelsStyle
    seqlevelsStyle(txdb) <- origTxDbLevelsStyle
  
    ## if the argument 'allTranscripts' is set to 'FALSE' then keep only once identical variants
    ## annotated with the same coding consequence to the same gene but in different tx
    if (!allTranscripts)
      GRanges_coding_uq <- GRanges_coding_uq[!duplicated(mcols(GRanges_coding_uq)[, c("QUERYID", "GENEID", "CONSEQUENCE")])]

    variantsVR_annotated_coding <- variantsVR_annotated_coding[GRanges_coding_uq$QUERYID]

    ## some TXID, CDSID and GENEID annotations are catched by
    ## VariantAnnotation::locateVariants() and not by VariantAnnotation::predictCoding()
    ## and the other way araound. In these cases we take the union of both annotations.
    ## when locateVariants() and predictCoding() disagree in the GENEID annotation then
    ## we take the one given by predictCoding() -SHOULD ASK ABOUT THIS IN THE DEVEL LIST-
    mask <- !is.na(variantsVR_annotated_coding$GENEID) & !is.na(GRanges_coding_uq$GENEID) &
            variantsVR_annotated_coding$GENEID != GRanges_coding_uq$GENEID

    variantsVR_annotated_coding$TXID[mask] <- GRanges_coding_uq$TXID[mask]
    variantsVR_annotated_coding$CDSID[mask] <- GRanges_coding_uq$CDSID[mask]
    variantsVR_annotated_coding$GENEID[mask] <- GRanges_coding_uq$GENEID[mask]

    mask <- is.na(variantsVR_annotated_coding$TXID)
    variantsVR_annotated_coding$TXID[mask] <- GRanges_coding_uq$TXID[mask]
    mask <- elementLengths(variantsVR_annotated_coding$CDSID) == 0
    variantsVR_annotated_coding$CDSID[mask] <- GRanges_coding_uq$CDSID[mask]
    mask <- is.na(variantsVR_annotated_coding$GENEID)
    variantsVR_annotated_coding$GENEID[mask] <- GRanges_coding_uq$GENEID[mask]

    ## add coding annotations from VariantAnnotation::predictCoding()
    variantsVR_annotated_coding$CDSLOC <- GRanges_coding_uq$CDSLOC
    variantsVR_annotated_coding$PROTEINLOC <- GRanges_coding_uq$PROTEINLOC
    variantsVR_annotated_coding$REFCODON <- GRanges_coding_uq$REFCODON
    variantsVR_annotated_coding$VARCODON <- GRanges_coding_uq$VARCODON
    variantsVR_annotated_coding$REFAA <- GRanges_coding_uq$REFAA
    variantsVR_annotated_coding$VARAA <- GRanges_coding_uq$VARAA
    variantsVR_annotated_coding$varAllele <- GRanges_coding_uq$varAllele
    variantsVR_annotated_coding$CONSEQUENCE <- GRanges_coding_uq$CONSEQUENCE
  
    ## annotate start and end positions of CDS (to determine start and stop gains and losses)
    cdsinfo <- select(txdb, keys=unique(as.character(unlist(variantsVR_annotated_coding$CDSID, use.names=FALSE))),
                      columns=c("CDSSTART", "CDSEND"), keytype="CDSID")
    mt <- match(as.character(unlist(variantsVR_annotated_coding$CDSID, use.names=FALSE)), cdsinfo$CDSID)
    variantsVR_annotated_coding$CDSSTART <- relist(cdsinfo$CDSSTART[mt], variantsVR_annotated_coding$CDSID)
    variantsVR_annotated_coding$CDSEND <- relist(cdsinfo$CDSEND[mt], variantsVR_annotated_coding$CDSID)

    ## annotate codon usage difference in synonymous mutations

    variantsVR_annotated_coding$CUREF <- NA_real_
    variantsVR_annotated_coding$CUALT <- NA_real_

    if (any(variantsVR_annotated_coding$CONSEQUENCE == "synonymous")){
    
      message("Annotating codon usage frequencies in coding synonymous variants")
    
      ## create a mask for those variants that are codying and synonymous

      mask <- variantsVR_annotated_coding$CONSEQUENCE %in% "synonymous"
    
      ## extract the reference codons and alt codons from the DNStringSet
   
      ref_codons <- unname(as.character(variantsVR_annotated_coding$REFCODON[mask]))
      alt_codons <- unname(as.character(variantsVR_annotated_coding$VARCODON[mask]))

      ## obtain the vector with the frequencies for the reference and alternative codons   
      ref_codons_numeric <- codonusageTable[ref_codons]
      alt_codons_numeric <- codonusageTable[alt_codons]
    
      ## push the numerical vector to its corresponding variant to compare on CUREF and CUALT columns
      variantsVR_annotated_coding$CUREF[mask] <- unname(ref_codons_numeric)
      variantsVR_annotated_coding$CUALT[mask] <- unname(alt_codons_numeric)
        
    }
  
  }

  ## consolidate the annotations on coding variants, obatined with 'VariantAnnotation::predictCoding()',
  ##  with the rest of non-coding variants in a single VRanges object, replacing the one named
  ## 'variantsVR_annotated'
  variantsVR_annotated_noncoding <- variantsVR_annotated[variantsVR_annotated$LOCATION != "coding"]
  n.noncoding <- length(variantsVR_annotated_noncoding)
  dummyDF <- DataFrame(CDSLOC=IRanges(start=rep(-1, n.noncoding), end=rep(-1, n.noncoding)),
                       PROTEINLOC=IntegerList(as.list(rep(NA, n.noncoding))),
                       REFCODON=DNAStringSet(rep("", n.noncoding)),
                       VARCODON=DNAStringSet(rep("", n.noncoding)),
                       REFAA=AAStringSet(rep("", n.noncoding)),
                       VARAA=AAStringSet(rep("", n.noncoding)),
                       varAllele=DNAStringSet(rep("", n.noncoding)),
                       CONSEQUENCE=factor(rep(NA, n.noncoding), levels=c("nonsynonymous", "synonymous")),
                       CDSSTART=IntegerList(as.list(rep(NA, n.noncoding))),
                       CDSEND=IntegerList(as.list(rep(NA, n.noncoding))),
                       CUREF=rep(NA_real_,n.noncoding),
                       CUALT=rep(NA_real_,n.noncoding))
  
  mcols(variantsVR_annotated_noncoding) <- cbind(mcols(variantsVR_annotated_noncoding), dummyDF)

  if (length(variantsVR_annotated_coding) > 0) {
    variantsVR_annotated <- sort(c(variantsVR_annotated_coding, variantsVR_annotated_noncoding))
  } else
    variantsVR_annotated <- sort(variantsVR_annotated_noncoding)

  ## annotate start and end positions of TX (to determine pre-mRNA positions)
  ## FIXME: with intergenic variants
  txinfo <- select(txdb, keys=unique(as.character(unlist(variantsVR_annotated$TXID, use.names=FALSE))),
                   columns=c("TXSTART", "TXEND"), keytype="TXID")
  mt <- match(as.character(unlist(variantsVR_annotated$TXID, use.names=FALSE)), txinfo$TXID)
  variantsVR_annotated$TXSTART <- txinfo$TXSTART[mt]
  variantsVR_annotated$TXEND <- txinfo$TXEND[mt]
  
  #############################################################
  ##                                                         ##
  ## ANNOTATE SPLICE SITES IN SYNONYMOUS & INTRONIC VARIANTS ##
  ##                                                         ##
  #############################################################
  
  ## add metadata columns in 'variantsVR_annotated' for cryptic ss annotations
  ## this should be optional once the shiny app is aware about present/absent annotations
  dummyDF <- DataFrame(SCORE5ssREF=rep(NA_real_, length(variantsVR_annotated)),
                       SCORE5ssALT=rep(NA_real_, length(variantsVR_annotated)),
                       SCORE5ssPOS=rep(NA_real_, length(variantsVR_annotated)),
                       SCORE3ssREF=rep(NA_real_, length(variantsVR_annotated)),
                       SCORE3ssALT=rep(NA_real_, length(variantsVR_annotated)),
                       SCORE3ssPOS=rep(NA_real_, length(variantsVR_annotated)))

  if (length(spliceSiteMatrices) == 2)
    dummyDF <- .scoreSpliceSiteVariants(variantsVR_annotated, spliceSiteMatrices, bsgenome, BPPARAM)
  mcols(variantsVR_annotated) <- cbind(mcols(variantsVR_annotated), dummyDF)

  ###############################
  ##                           ##
  ##  GENE-CENTRIC ANNOTATIONS ##
  ##                           ##
  ###############################

  mcols(variantsVR_annotated) <- cbind(mcols(variantsVR_annotated),
                                       annotateVariants(orgdb, variantsVR_annotated, param))

  #####################################
  ##                                 ##
  ##  TRANSCRIPT-CENTRIC ANNOTATIONS ##
  ##                                 ##
  #####################################

  mcols(variantsVR_annotated) <- cbind(mcols(variantsVR_annotated),
                                       annotateVariants(txdb, variantsVR_annotated, param))

  ########################
  ##                    ##
  ##    DESCRIPTION     ##
  ##                    ##
  ########################

  mcols(variantsVR_annotated) <- cbind(mcols(variantsVR_annotated), variantHGVS(variantsVR_annotated))
                                      
  ########################
  ##                    ##
  ##     AminoAcid      ##
  ##                    ##
  ########################

  mcols(variantsVR_annotated) <- cbind(mcols(variantsVR_annotated),
                                       aminoAcidChanges(variantsVR_annotated, radicalAAchangeMatrix))

  ########################
  ##                    ##
  ##  OTHER ANNOTATIONS ##
  ##                    ##
  ########################
  
  for (i in seq(along=otherAnnotations)) {
    message(sprintf("Annotating with %s", names(otherAnnotations)[i]))
    mcols(variantsVR_annotated) <- cbind(mcols(variantsVR_annotated),
                                         annotateVariants(otherAnnotations[[i]],
                                                          variantsVR_annotated,
                                                          param,
                                                          BPPARAM=BPPARAM))
  }

  ## this seems to be nicely transferred during the coercion to VRanges
  ## restore the seqinfo data
  ## seqinfo(variantsVR_annotated) <- seqinfo(variantsVR)

  return(variantsVR_annotated)
}



###########
## Annotate dbSNP identifiers
####

setMethod("annotateVariants", signature(annObj="SNPlocs"),
          function(annObj, variantsVR, param, BPPARAM=bpparam("SerialParam")) {
            if (!"TYPE" %in% colnames(mcols(variantsVR))) {
              stop("Variant type (SNV, Insertion, Deletion, MNV, Delins) has not been annotated.")
            }
            seqlevelsStyle(variantsVR) <- seqlevelsStyle(annObj)
            masksnp <- variantsVR$TYPE == "SNV"
            rsids <- rep(NA_character_, times=length(variantsVR))
            if (any(masksnp)) {
              rsids_list <- .loc2SNPid(annObj, variantsVR[masksnp], BPPARAM=BPPARAM)
              elen <- elementLengths(rsids_list)
              rsids[masksnp][elen == 1] <- as.character(rsids_list[elen == 1])

              ## paste together multiple dbSNP identifiers
              rsids[masksnp][elen > 1] <- unlist(bplapply(rsids_list[elen > 1], paste, collapse=", ", BPPARAM=BPPARAM),
                                               use.names=FALSE)
              rsids[rsids == ""] <- NA_character_
            }
            return(DataFrame(dbSNP=rsids))
          })

setMethod("annotateVariants", signature(annObj="XtraSNPlocs"),
          function(annObj, variantsVR, param, BPPARAM=bpparam("SerialParam")) {
            if (!"TYPE" %in% colnames(mcols(variantsVR))) {
              stop("Variant type (SNV, Insertion, Deletion, MNV, Delins) has not been annotated.")
            }
            seqlevelsStyle(variantsVR) <- seqlevelsStyle(annObj)
            maskxtrasnp <- variantsVR$TYPE != "SNV"
            rsids <- rep(NA_character_, times=length(variantsVR))
            if (any(maskxtrasnp)) {
              rsids_list <- .loc2XtraSNPid(annObj, variantsVR[maskxtrasnp], BPPARAM=BPPARAM)
              elen <- elementLengths(rsids_list)
              rsids[maskxtrasnp][elen == 1] <- as.character(rsids_list[elen == 1])

              ## paste together multiple dbSNP identifiers
              rsids[maskxtrasnp][elen > 1] <- unlist(bplapply(rsids_list[elen > 1], paste, collapse=", ", BPPARAM=BPPARAM),
                                                     use.names=FALSE)
              rsids[rsids == ""] <- NA_character_
            }
            return(DataFrame(dbSNP=rsids))
          })

###########
## Annotate PolyPhen2 predictions
####

setMethod("annotateVariants", signature(annObj="PolyPhenDb"),
          function(annObj, variantsVR, param, coding=TRUE, BPPARAM=bpparam("SerialParam")) {
            PolyPhen2 <- rep(NA_character_, length(variantsVR))
            if (!coding)
              return(DataFrame(PolyPhen2=PolyPhen2))

            rsids <- unique(variantsVR$dbSNP)
            rsids <- rsids[!is.na(rsids)]
            if (length(rsids) > 0) {
              pp <- select(annObj, keys=rsids,
                           cols=c("TRAININGSET", "PREDICTION", "PPH2PROB"))
              mask <- pp$TRAININGSET == "humvar"
              pphumvar <- pp[mask, ]
              mt <- match(variantsVR$dbSNP, pphumvar$RSID)
              PolyPhen2[!is.na(mt)] <- pphumvar$PREDICTION[mt[!is.na(mt)]]
            }
            DataFrame(PolyPhen2=PolyPhen2)
          })

############
## Annotate PROVEAN predictions (former SIFT)
#####

## USE VARID WHERE dbSNP IS MISSING !!!
setMethod("annotateVariants", signature(annObj="PROVEANDb"),
          function(annObj, variantsVR, param, coding=TRUE, BPPARAM=bpparam("SerialParam")) {
            PROVEAN <- rep(NA_character_, length(variantsVR))
            if (!coding)
              return(DataFrame(PROVEAN=PROVEAN))

            rsids <- unique(variantsVR$dbSNP)
            rsids <- rsids[!is.na(rsids)]
            if (length(rsids) > 0) {
              rsids <- gsub("rs", "", rsids)
              pv <- select(annObj, keys=rsids, columns="PROVEANPRED")
              mt <- match(gsub("rs", "", variantsVR$dbSNP), pv$DBSNPID)
              PROVEAN[!is.na(mt)] <- pv$PROVEANPRED[mt[!is.na(mt)]]
            }
            DataFrame(PROVEAN=PROVEAN)
          })

############
## Annotate MAF values
#####

setMethod("annotateVariants", signature(annObj="MafDb"),
          function(annObj, variantsVR, param, BPPARAM=bpparam("SerialParam")) {

            ## get the MAF columns
            mafCols <- knownVariantsMAFcols(annObj) ## assumes all MAF column names contain 'AF'

            mafValues <- matrix(NA, nrow=length(variantsVR), ncol=length(mafCols),
                                dimnames=list(NULL, mafCols))
            varIDs <- variantsVR$dbSNP     ## fetch first by the annotated dbSNP identifier
            uniqVarIDs <- unique(varIDs)
            uniqVarIDs <- uniqVarIDs[!is.na(uniqVarIDs)]
            if (length(varIDs) > 0) {
              uniqMAFvalues <- snpid2maf(annObj, uniqVarIDs)
              mt <- match(varIDs, uniqMAFvalues$varID)
              mafValues[!is.na(mt), ] <- as.matrix(uniqMAFvalues[mt[!is.na(mt)], mafCols, drop=FALSE])

              ## for missing entries then fetch by given identifier
              varIDs <- variantsVR$VARID
              if (!is.null(varIDs) && any(!is.na(varIDs))) {
                maskNotFound <- apply(uniqMAFvalues[, mafCols], 1, function(x) all(is.na(x)))
                missingdbsnpIDs <- unique(c(varIDs[is.na(variantsVR$dbSNP)],
                                            varIDs[!is.na(variantsVR$dbSNP) &
                                                   variantsVR$dbSNP %in% uniqMAFvalues$varID[maskNotFound]]))
                missingdbsnpIDs <- missingdbsnpIDs[!is.na(missingdbsnpIDs)]
                if (length(missingdbsnpIDs) > 0) {
                  uniqMAFvalues <- snpid2maf(annObj, missingdbsnpIDs)
                  mt <- match(varIDs, uniqMAFvalues$varID)
                  mafValues[!is.na(mt), ] <- as.matrix(uniqMAFvalues[mt[!is.na(mt)], mafCols, drop=FALSE])
                }
              }
            }

            colnames(mafValues) <- paste0(colnames(mafValues), annObj$tag) ## tag MAF columns with their data source

            DataFrame(mafValues)
          })

############
## Annotate organism-level gene-centric features
#####

###########
## Get HGNC gene symbol
####

setMethod("annotateVariants", signature(annObj="OrgDb"),
          function(annObj, variantsVR, param, BPPARAM=bpparam("SerialParam")) {

            genelevel_annot <- DataFrame(GENE=character(), OMIM=character())
            geneIDs <- variantsVR$GENEID
            geneKeytype <- param$geneKeytype
            maskNAs <- is.na(geneIDs)
            if (length(geneIDs) > 0) {
              ## if input IDs are NAs output should also be NAs and avoid querying malformed keys
              genelevel_annot <- DataFrame(GENE=rep(NA_character_, times=length(geneIDs)),
                                           OMIM=rep(NA_character_, times=length(geneIDs)))
              if (sum(!maskNAs) > 0) {
                if (is.na(geneKeytype)) {
                  geneKeytype <- "ENTREZID"
                  if (substr(geneIDs[!maskNAs][1], 1, 4) == "ENSG")
                    geneKeytype <- "ENSEMBL"
                  else if (substr(geneIDs[!maskNAs][1], 1, 3) %in% c("NM_", "NP_", "NR_", "XM_", "XP", "XR_", "YP_"))
                    geneKeytype <- "REFSEQ"
                }
                uniqEntrezIDs <- unique(geneIDs[!maskNAs])
                res <- select(annObj, keys=as.character(uniqEntrezIDs), columns=c("SYMBOL", "OMIM"), keytype=geneKeytype)
                symxgeneID <- sapply(split(res$SYMBOL, res[[geneKeytype]]),
                                     function(x) paste(unique(x), collapse=", "))
                omimxgeneID <- sapply(split(res$OMIM, res[[geneKeytype]]),
                                      function(x) {
                                        if (all(!is.na(x))) x <- paste(unique(x), collapse=", ") ; x
                                      })
                genelevel_annot[!maskNAs, ] <- DataFrame(GENE=symxgeneID[geneIDs[!maskNAs]],
                                                         OMIM=omimxgeneID[geneIDs[!maskNAs]])
              }
            }
            genelevel_annot
          })

setMethod("annotateVariants", signature(annObj="TxDb"),
          function(annObj, variantsVR, param, BPPARAM=bpparam("SerialParam")) {

            txlevel_annot <- DataFrame(TXNAME=character())
            txIDs <- variantsVR$TXID
            maskNAs <- is.na(txIDs)
            if (length(txIDs) > 0) {
              ## if input IDs are NAs output should also be NAs and avoid querying malformed keys
              txlevel_annot <- DataFrame(TXNAME=rep(NA_character_, times=length(txIDs)))
              if (sum(!maskNAs) > 0) {
                uniqTxIDs <- unique(txIDs[!maskNAs])
                res <- select(annObj, keys=as.character(uniqTxIDs), columns="TXNAME", keytype="TXID")
                txnamextxID <- sapply(split(res$TXNAME, res$TXID),
                                       function(x) paste(unique(x), collapse=", "))
                txlevel_annot[!maskNAs, ] <- DataFrame(TXNAME=txnamextxID[txIDs[!maskNAs]])
              }
            }
            txlevel_annot
          })

############
## Annotate phastCons values
#####

setMethod("annotateVariants", signature(annObj="PhastConsDb"),
          function(annObj, variantsVR, param, BPPARAM=bpparam("SerialParam")) {

            sco <- scores(annObj, variantsVR)

            DataFrame(phastCons=sco)
          })

############
## Annotate gene phylostratum
#####

setMethod("annotateVariants", signature(annObj="GenePhylostrataDb"),
          function(annObj, variantsVR, param, BPPARAM=bpparam("SerialParam")) {

            gps <- genePhylostratum(annObj, variantsVR$GENEID)

            DataFrame(GenePhylostratumTaxID=gps$TaxID,
                      GenePhylostratumIndex=gps$OldestPhylostratum,
                      GenePhylostratum=gps$Description)
          })


###########
## Read matrix of radical versus conservative amino acid substitutions
#####

readAAradicalChangeMatrix <- function(file) {

  ## read amino acid properties and build matrix of radical changes
  aaProperties <- read.delim(file, comment.char="#")
  aaCodes <- rownames(aaProperties)

  if (is.null(aaCodes) || any(is.na(match(aaCodes, c(LETTERS, "*")))))
    stop(sprintf("Some row names in %s do not correspond to amino acid IUPAC one-letter codes.", file))

  classColumns <- apply(aaProperties[, 3:ncol(aaProperties)], 2, class)
  if (any(classColumns != "logical"))
    stop(sprintf("Column(s) %s in %s does not form a logical vector or TRUE/FALSE values.",
                 paste(which(classColumns != "logical"), collapse=", "), file))

  radicalAAchanges <- matrix(TRUE, nrow=nrow(aaProperties), ncol=nrow(aaProperties),
                             dimnames=list(aaCodes, aaCodes))

  ## amino acid changes within a property group are *not* radical, also called conservative.
  for (i in 3:ncol(aaProperties)) { ## properties start at the third column
    prop <- aaProperties[[i]]
    conservativeChanges <- t(combn(aaCodes[prop], 2))
    radicalAAchanges[conservativeChanges] <- FALSE
  }
  radicalAAchanges <- radicalAAchanges & t(radicalAAchanges) ## made it symmetric
  diag(radicalAAchanges) <- FALSE ## the change of an amino acid by itself is not radical

  ## consider the change condition with respect to an unknown amino acid as missing
  radicalAAchanges["X", ] <- radicalAAchanges[, "X"] <- NA

  radicalAAchanges
}



###
### PRIVATE FUNCTIONS
###

typeOfVariants <- function(variantsVR) {

  type <- factor(levels=c("SNV", "Insertion", "Deletion", "MNV", "Delins"))
  if (length(variantsVR) > 0) {
    type <- factor(rep("SNV", times=length(variantsVR)),
                   levels=c("SNV", "Insertion", "Deletion", "MNV", "Delins"))
    type[isInsertion(variantsVR)] <- "Insertion"
    type[isDeletion(variantsVR)] <- "Deletion"
    type[isSubstitution(variantsVR) & !isSNV(variantsVR)] <- "MNV"
    type[isDelins(variantsVR)] <- "Delins"
  }

  DataFrame(TYPE=type)
}

## this function tries to provide a variant description following
## HGVS nomenclature at http://www.hgvs.org/mutnomen
## this is a very first version covering only coding descriptions
## hopefully it will become more comprehensive in the near future
variantHGVS <- function(variantsVR) {

  if (length(variantsVR) == 0)
    return(DataFrame(HGVSg=character(), HGVSc=character(), HGVSp=character()))

  if (!all(c("LOCATION", "TYPE", "cDNALOC") %in% colnames(mcols(variantsVR))))
    stop("Metadata columns LOCATION, TYPE and cDNALOC should be annotated before calling variantDescription().")

  pDesc <- gDesc <- cDesc <- refAllele <- altAllele <- rep(NA_character_, length(variantsVR))
  locAllele <- rep(NA_integer_, length(variantsVR))

  ## HGVS coding annotations

  maskCoding <- variantsVR$LOCATION == "coding"

  if (any(maskCoding))
    ## for coding variants we use the 'varAllele' column which is already adjusted for strand
    refAllele[maskCoding] <- .adjustForStrandSense(variantsVR[maskCoding], ref(variantsVR)[maskCoding])

    ## for non-coding variants we have to adjust for strand both, reference and alternative alleles
    ## THIS IS PROBABLY REDUNDANT AS THE VRanges CONSTRUCTOR ALREADY ADJUSTS FOR THIS (???)
  refAllele[!maskCoding] <- .adjustForStrandSense(variantsVR[!maskCoding], ref(variantsVR)[!maskCoding])
  altAllele <- as.character(variantsVR$varAllele)
  altAllele[!maskCoding] <- .adjustForStrandSense(variantsVR[!maskCoding], alt(variantsVR)[!maskCoding])

  locStartAllele <- as.integer(start(variantsVR$CDSLOC))
  locEndAllele <- as.integer(end(variantsVR$CDSLOC))
  widthAllele <- as.integer(width(variantsVR$CDSLOC))

  ## SNVs
  mask <- maskCoding & variantsVR$TYPE == "SNV"
  cDesc[mask] <- sprintf("c.%d%s>%s", locStartAllele[mask], refAllele[mask], altAllele[mask])

  ## Insertions
  mask <- maskCoding & variantsVR$TYPE == "Insertion"
  cDesc[mask] <- sprintf("c.%d_%dins%s", locStartAllele[mask], locStartAllele[mask]+1, altAllele[mask])

  ## Deletions
  mask <- maskCoding & variantsVR$TYPE == "Deletion" & widthAllele == 1
  cDesc[mask] <- sprintf("c.%ddel%s", locStartAllele[mask], altAllele[mask])

  mask <- maskCoding & variantsVR$TYPE == "Deletion" & widthAllele > 1
  cDesc[mask] <- sprintf("c.%d_%ddel%s", locStartAllele[mask], locEndAllele[mask], altAllele[mask])

  ## Deletions-insertions
  mask <- maskCoding & variantsVR$TYPE == "Delins"
  cDesc[mask] <- sprintf("c.%d_%ddelins%s", locStartAllele[mask], locEndAllele[mask], altAllele[mask])

  ## HGVS genomic annotations

  locStartAllele <- locEndAllele <- rep(NA_integer_, times=length(variantsVR))

  mask <- variantsVR$LOCATION != "intergenic"
  if (any(mask)) {
    locStartAllele[mask] <- ifelse(strand(variantsVR)[mask] == "+", as.integer(start(variantsVR)[mask]) - variantsVR$TXSTART[mask] + 1L,
                                   variantsVR$TXEND[mask] - as.integer(end(variantsVR)[mask]) + 1L)
    locEndAllele[mask] <- ifelse(strand(variantsVR)[mask] == "+", as.integer(end(variantsVR)[mask]) - variantsVR$TXSTART[mask] + 1L,
                                 variantsVR$TXEND[mask] - as.integer(start(variantsVR)[mask]) + 1L)
  }

  if (any(!mask)) { ## for intergenic variants just use their position as given
    locStartAllele[!mask] <- as.integer(start(variantsVR)[!mask])
    locEndAllele[!mask] <- as.integer(end(variantsVR)[!mask])
  }
  widthAllele <- as.integer(width(variantsVR))
  
  ## SNVs
  mask <- variantsVR$TYPE == "SNV"
  gDesc[mask] <- sprintf("g.%d%s>%s", locStartAllele[mask], refAllele[mask], altAllele[mask])

  ## Insertions
  mask <- variantsVR$TYPE == "Insertion"
  gDesc[mask] <- sprintf("g.%d_%dins%s", locStartAllele[mask], locEndAllele[mask], altAllele[mask])

  ## Deletions
  mask <- variantsVR$TYPE == "Deletion" & widthAllele == 1
  gDesc[mask] <- sprintf("g.%ddel%s", locStartAllele[mask], altAllele[mask])
  mask <- variantsVR$TYPE == "Deletion" & widthAllele > 1
  gDesc[mask] <- sprintf("g.%d_%ddel%s", locStartAllele[mask], locEndAllele[mask], altAllele[mask])

  ## Deletions-insertions
  mask <- variantsVR$TYPE == "Delins"
  gDesc[mask] <- sprintf("g.%d_%ddelins%s", locStartAllele[mask], locEndAllele[mask], altAllele[mask])

  ## HGVS protein annotations

  locStartAllele <- as.integer(sapply(variantsVR$PROTEINLOC, "[", 1))
  locEndAllele <- as.integer(sapply(variantsVR$PROTEINLOC, "[", 2))
  mask <- !is.na(locStartAllele) & is.na(locEndAllele)
  locEndAllele[mask] <- locStartAllele[mask]
  widthAllele <- locEndAllele - locStartAllele + 1L

  ## SNVs
  mask <- maskCoding & variantsVR$TYPE == "SNV"
  pDesc[mask] <- sprintf("p.%d%s>%s", locStartAllele[mask], variantsVR$REFAA[mask], variantsVR$VARAA[mask])

  ## Insertions
  mask <- maskCoding & variantsVR$TYPE == "Insertion"
  pDesc[mask] <- sprintf("p.%d_%dins%s", locStartAllele[mask], locStartAllele[mask] + 1L, variantsVR$VARAA[mask])

  ## Deletions
  mask <- maskCoding & variantsVR$TYPE == "Deletion" & widthAllele == 1
  pDesc[mask] <- sprintf("p.%ddel%s", locStartAllele[mask], variantsVR$VARAA[mask])

  mask <- maskCoding & variantsVR$TYPE == "Deletion" & widthAllele > 1
  pDesc[mask] <- sprintf("p.%d_%ddel%s", locStartAllele[mask], locEndAllele[mask], variantsVR$VARAA[mask])

  ## Deletions-insertions
  mask <- maskCoding & variantsVR$TYPE == "Delins"
  pDesc[mask] <- sprintf("p.%d_%ddelins%s", locStartAllele[mask], locEndAllele[mask], variantsVR$VARAA[mask])
  
  DataFrame(HGVSg=gDesc, HGVSc=cDesc, HGVSp=pDesc)
}

## adapted from http://permalink.gmane.org/gmane.science.biology.informatics.conductor/48456
.loc2SNPid <- function(SNPlocsObj, locs, BPPARAM=bpparam("SerialParam")) {

  if (!is(locs, "GRanges"))
    stop("'locs' must be a GRanges object")
  if (!all(width(locs) == 1L))
    stop("all ranges in 'locs' must be of width 1")
  common_seqlevels <- intersect(seqlevels(locs), names(snpcount(SNPlocsObj)))
  if (length(common_seqlevels) == 0L)
    stop("chromosome names (a.k.a. seqlevels) in 'locs' don't seem to ",
          "be\n  compatible with the chromosome names in the SNPlocs ",
          "package. Maybe they\n  use a different naming convention? ",
          "If that's the case then you first need\n  to rename the ",
          "seqlevels in 'locs'. See '?seqlevels' for how to do this.")
  f <- as.factor(seqnames(locs))
  locs_by_chrom <- split(start(locs), f)
  rsids_by_chrom <- bplapply(seq_along(locs_by_chrom),
                             function(i) {
                               seqname <- levels(f)[i]
                               ## locs2 <- locs_by_chrom[[i]]
                               nlocs_by_chrom <- length(locs_by_chrom[[i]])
                               locs2 <- GRanges(Rle(seqname, nlocs_by_chrom),
                                                IRanges(start=locs_by_chrom[[i]],
                                                        width=rep(1, times=nlocs_by_chrom)))
                               ans2 <- vector("list", length=length(locs2))
                               if (length(locs2) == 0L || !(seqname %in% common_seqlevels))
                                   return(ans2)
                               locs3 <- snplocs(SNPlocsObj, seqname, as.GRanges=TRUE)
                               hits <- findOverlaps(locs2, locs3) ## findOverlaps on a GRanges faster than findMatches on a vector
                               ## hits <- findMatches(locs2, locs3$loc)
                               if (length(hits) > 0) {
                                   rsids <- paste0("rs", locs3$RefSNP_id[subjectHits(hits)])
                                   q_hits <- queryHits(hits)
                                   tmp <- split(rsids, q_hits)
                                   ans2[as.integer(names(tmp))] <- tmp
                               }
                               ans2
                           }, BPPARAM=BPPARAM)
  CharacterList(unsplit(rsids_by_chrom, f))
}

## adapted from http://permalink.gmane.org/gmane.science.biology.informatics.conductor/48456
.loc2XtraSNPid <- function(XtraSNPlocsObj, locs, BPPARAM=bpparam("SerialParam")) {

  if (!is(locs, "GRanges"))
    stop("'locs' must be a GRanges object")

  mcols(locs) <- NULL
  locs <- as(locs, "GRanges") ## when 'locs' is a 'VRanges' object

  common_seqlevels <- intersect(seqlevels(locs), names(snpcount(XtraSNPlocsObj)))
  if (length(common_seqlevels) == 0L)
    stop("chromosome names (a.k.a. seqlevels) in 'locs' don't seem to ",
          "be\n  compatible with the chromosome names in the SNPlocs ",
          "package. Maybe they\n  use a different naming convention? ",
          "If that's the case then you first need\n  to rename the ",
          "seqlevels in 'locs'. See '?seqlevels' for how to do this.")
  f <- as.factor(seqnames(locs))
  locs_by_chrom <- split(locs, f)
  rsids_by_chrom <- bplapply(locs_by_chrom,
                             function(locs2) {
                               ans2 <- vector("list", length=length(locs2))
                               if (length(locs2) == 0L)
                                   return(ans2)

                               seqname <- as.character(seqnames(locs2)[1])
                               if (!seqname %in% common_seqlevels)
                                   return(ans2)
                               
                               locs3 <- snpsByOverlaps(XtraSNPlocsObj, locs2, type="any",
                                                       columns="RefSNP_id")
                               hits <- findOverlaps(locs2, locs3)

                               if (length(hits) > 0) {
                                   rsids <- locs3$RefSNP_id[subjectHits(hits)]
                                   q_hits <- queryHits(hits)
                                   tmp <- split(rsids, q_hits)
                                   ans2[as.integer(names(tmp))] <- tmp
                               }
                               ans2
                           }, BPPARAM=BPPARAM)
  CharacterList(unsplit(rsids_by_chrom, f))
}


## annotate the amino acid substitution in a XposX format with X being the
## IUPAC amino acid code and pos the position along the protein, and whether
## the change in amino acid can be considered chemically radical or conservative
## frameshifts and nonsense changes are considered directly radical
aminoAcidChanges <- function(variantsVR, rAAch) {
  aachange <- aachangetype <- rep(NA_character_, length(variantsVR))

  ## annotate non-synonymous coding changes
  whnonsyn <- which(variantsVR$CONSEQUENCE == "nonsynonymous")
  locaa <- rep(NA_character_, length(length(whnonsyn)))
  refaa <- as.character(variantsVR[whnonsyn]$REFAA)
  altaa <- as.character(variantsVR[whnonsyn]$VARAA)
  elen <- elementLengths(variantsVR[whnonsyn]$PROTEINLOC)
  
  locaa[elen == 1] <- as.character(unlist(variantsVR[whnonsyn]$PROTEINLOC[elen == 1], use.names=FALSE))

  ## location of a multiple amino acid replacement is denoted by its position range
  locaa[elen > 1] <- sapply(variantsVR[whnonsyn]$PROTEINLOC[elen > 1],
                                      function(x) paste(range(x), collapse="-"))

  ## denote amino acid changes by concatenating reference amino acid, position
  ## and alternative amino acid
  aachange[whnonsyn] <- paste0(refaa, locaa, altaa)

  ## on single amino acid replacements evaluate the type of change (radical or conservative)
  ## using the input logical matrix 'rAAch' whose cells should be set to true when the
  ## change is radical
  masksinglechanges <- nchar(refaa) == 1 & nchar(altaa) == 1
  aachangetype[whnonsyn][masksinglechanges] <- ifelse(rAAch[cbind(refaa[masksinglechanges], altaa[masksinglechanges])],
                                                      "Radical", "Conservative")

  ## annotate synonymous coding no-changes as "Conservative" and
  ## to indicate their change position along the protein sequence
  whsyn <- which(variantsVR$CONSEQUENCE == "synonymous")
  locaa <- rep(NA_character_, length(length(whsyn)))
  refaa <- as.character(variantsVR[whsyn]$REFAA)
  altaa <- as.character(variantsVR[whsyn]$VARAA)
  locaa <- as.character(unlist(variantsVR[whsyn]$PROTEINLOC, use.names=FALSE))
  aachange[whsyn] <- paste0(refaa, locaa, altaa)
  aachangetype[whsyn] <- "Conservative"

  ## annotate frameshift changes as "Radical"
  whframeshift <- which(variantsVR$CONSEQUENCE == "frameshift")
  locaa <- rep(NA_character_, length(length(whframeshift)))
  refaa <- as.character(variantsVR[whframeshift]$REFAA)
  altaa <- as.character(variantsVR[whframeshift]$VARAA)
  elen <- elementLengths(variantsVR[whframeshift]$PROTEINLOC)
  locaa[elen == 1] <- as.character(unlist(variantsVR[whframeshift]$PROTEINLOC[elen == 1], use.names=FALSE))
  ## location of a multiple amino acid replacement is denoted by its position range
  locaa[elen > 1] <- sapply(variantsVR[whframeshift]$PROTEINLOC[elen > 1],
                                      function(x) paste(range(x), collapse="-"))
  aachange[whframeshift] <- paste0(refaa, locaa, altaa)
  aachangetype[whframeshift] <- "Radical"

  ## annotate nonsense changes as "Radical"
  whnonsense <- which(variantsVR$CONSEQUENCE == "nonsense")
  locaa <- rep(NA_character_, length(length(whnonsense)))
  refaa <- as.character(variantsVR[whnonsense]$REFAA)
  altaa <- as.character(variantsVR[whnonsense]$VARAA)
  locaa <- as.character(unlist(variantsVR[whnonsense]$PROTEINLOC, use.names=FALSE))
  aachange[whnonsense] <- paste0(refaa, locaa, altaa)
  aachangetype[whnonsense] <- "Radical"

  DataFrame(AAchange=aachange, AAchangeType=aachangetype)
}

.emptyAnnotations <- function() {
  DataFrame(LOCATION=factor(),
            LOCSTART=integer(),
            LOCEND=integer(),
            TXID=integer(),
            CDSID=integer(),
            GENEID=character(),
            REF=DNAStringSet(),
            ALT=DNAStringSetList(),
            TYPE=character(),
            dbSNP=character(),
            cDNALOC=integer(),
            varAllele=DNAStringSet(),
            CDSLOC=IRanges(),
            PROTEINLOC=IntegerList(),
            CONSEQUENCE=factor(),
            REFCODON=DNAStringSet(),
            VARCODON=DNAStringSet(),
            REFAA=AAStringSet(),
            VARAA=AAStringSet(),
            CUREF=numeric(),
            CUALT=numeric(),
            TXSTART=integer(),
            TXEND=integer(),
            SCORE5ssREF=numeric(),
            SCORE5ssALT=numeric(),
            SCORE5ssPOS=numeric(),
            SCORE3ssREF=numeric(),
            SCORE3ssALT=numeric(),
            SCORE3ssPOS=numeric(),
            GENE=character(),
            OMIM=character(),
            TXNAME=character(),
            HGVSg=character(),
            HGVSc=character(),
            HGVSp=character(),
            AAchange=character(),
            AAchangeType=character())
}

## function to score splice sites including variants. it produces scores for the splice site with
## the reference and alternative alleles. it assumes that the input variantsVR is a VRanges object
.scoreSpliceSiteVariants <- function(variantsVR, spliceSiteMatrices, bsgenome, BPPARAM=bpparam("SerialParam")) {

  ## adapt to sequence style and genome version from the input
  ## BSgenome object, thus assuming positions are based on the same
  ## genome even though might be differently specified (i.e., hg38 vs GRCh38)
  seqlevelsStyle(variantsVR) <- seqlevelsStyle(bsgenome)
  commonChr <- intersect(seqlevels(variantsVR), seqlevels(bsgenome))
  if (any(is.na(genome(variantsVR)))) {
    warning(sprintf("Assuming the genome build of the input variants is %s.", unique(genome(bsgenome)[commonChr])))
    genome(variantsVR) <- genome(bsgenome)
  } else if (any(genome(variantsVR)[commonChr] != genome(bsgenome)[commonChr])) {
    warning(sprintf("Assumming %s represent the same genome build.",
                    paste(c(unique(genome(variantsVR)[commonChr]), unique(genome(bsgenome)[commonChr])),
                          collapse=" and ")))
    genome(variantsVR) <- genome(bsgenome)
  }

  ## add metadata columns in 'variantsVR' for splice site score annotations
  dummyDF <- DataFrame(SCORE5ssREF=rep(NA_real_, length(variantsVR)),
                       SCORE5ssALT=rep(NA_real_, length(variantsVR)),
                       SCORE5ssPOS=rep(NA_integer_, length(variantsVR)),
                       SCORE3ssREF=rep(NA_real_, length(variantsVR)),
                       SCORE3ssALT=rep(NA_real_, length(variantsVR)),
                       SCORE3ssPOS=rep(NA_integer_, length(variantsVR)))

  wmDonorSites <- spliceSiteMatrices$wmDonorSites
  wmAcceptorSites <- spliceSiteMatrices$wmAcceptorSites

  ## annotated splice sites

  message("Scoring annotated 5' splice sites")

  ssSNVmask <- variantsVR$TYPE == "SNV" & variantsVR$LOCATION == "fiveSpliceSite"
  if (any(ssSNVmask)) {
    GRanges_annotSS <- variantsVR[ssSNVmask]

    wregion <- GRanges_annotSS$LOCEND[1] - GRanges_annotSS$LOCSTART[1] + 1
    if (wregion == width(wmDonorSites)) {

      # get alternative allele adjusted by strand, we need a DNAStringSetList to use replaceAt() below
      altAlleleStrandAdjusted <- DNAStringSetList(strsplit(alt(GRanges_annotSS), split="", fixed=TRUE))
      altAlleleStrandAdjusted <- .adjustForStrandSense(GRanges_annotSS, altAlleleStrandAdjusted)

      GRanges_annotSS <- GRanges(seqnames=seqnames(GRanges_annotSS),
                                 ranges=IRanges(GRanges_annotSS$LOCSTART, GRanges_annotSS$LOCEND),
                                 strand=strand(GRanges_annotSS),
                                 POS=start(GRanges_annotSS) - GRanges_annotSS$LOCSTART + 1)

      # retrieve region of the splice site including the reference allele
      GRanges_annotSS_REF_strings <- getSeq(bsgenome, GRanges_annotSS)

      # replace the variant by the alternate allele
      GRanges_annotSS_ALT_strings <- replaceAt(GRanges_annotSS_REF_strings,
                                               IRangesList(start=as.list(GRanges_annotSS$POS),
                                                           end=as.list(GRanges_annotSS$POS)),
                                               altAlleleStrandAdjusted)

      # score REF alleles for donor splice sites
      GRanges_annotSS_REF_scores <- bpvec(X=GRanges_annotSS_REF_strings,
                                          FUN=wmScore, object=wmDonorSites, BPPARAM=BPPARAM)

      # score ALT alleles for donor splice sites
      GRanges_annotSS_ALT_scores <- bpvec(X=GRanges_annotSS_ALT_strings,
                                          FUN=wmScore, object=wmDonorSites, BPPARAM=BPPARAM)

      # store the position of the allele respect to the position of the dinucleotide GT whose nucleotides occur at pos 1 and 2
      relposdonor <- seq(-conservedPositions(wmDonorSites)[1]+1, -1, by=1)
      relposdonor <- c(relposdonor, seq(1, width(wmDonorSites)-length(relposdonor)))

      SCOREss <- DataFrame(SCORE5ssREF=round(GRanges_annotSS_REF_scores, digits=2),
                           SCORE5ssALT=round(GRanges_annotSS_ALT_scores, digits=2),
                           SCORE5ssPOS=relposdonor[GRanges_annotSS$POS],
                           SCORE3ssREF=rep(NA_real_, length(GRanges_annotSS)),
                           SCORE3ssALT=rep(NA_real_, length(GRanges_annotSS)),
                           SCORE3ssPOS=rep(NA_integer_, length(GRanges_annotSS)))

      ## incorporate the cryptic splice site annotations on synonymous variants into 'variantsVR'
      dummyDF[ssSNVmask, ] <- SCOREss
    } else
      warning(sprintf("Width of the 5' splice site region (%d) is not equal to the width of the weight matrix for donor sites (%d).",
                      wregion, width(wmDonorSites)))
  }

  message("Scoring annotated 3' splice sites")

  ssSNVmask <- variantsVR$TYPE == "SNV" & variantsVR$LOCATION == "threeSpliceSite"
  if (any(ssSNVmask)) {
    GRanges_annotSS <- variantsVR[ssSNVmask]

    wregion <- GRanges_annotSS$LOCEND[1] - GRanges_annotSS$LOCSTART[1] + 1
    if (wregion == width(wmAcceptorSites)) {

      # get alternative allele adjusted by strand, we need a DNAStringSetList to use replaceAt() below
      altAlleleStrandAdjusted <- DNAStringSetList(strsplit(alt(GRanges_annotSS), split="", fixed=TRUE))
      altAlleleStrandAdjusted <- .adjustForStrandSense(GRanges_annotSS, altAlleleStrandAdjusted)

      GRanges_annotSS <- GRanges(seqnames=seqnames(GRanges_annotSS),
                                 ranges=IRanges(GRanges_annotSS$LOCSTART, GRanges_annotSS$LOCEND),
                                 strand=strand(GRanges_annotSS),
                                 POS=start(GRanges_annotSS) - GRanges_annotSS$LOCSTART + 1)

      # retrieve region of the splice site including the reference allele
      GRanges_annotSS_REF_strings <- getSeq(bsgenome, GRanges_annotSS)

      # replace the variant by the alternate allele
      GRanges_annotSS_ALT_strings <- replaceAt(GRanges_annotSS_REF_strings,
                                               IRangesList(start=as.list(GRanges_annotSS$POS),
                                                           end=as.list(GRanges_annotSS$POS)),
                                               altAlleleStrandAdjusted)

      # score REF alleles for donor splice sites
      GRanges_annotSS_REF_scores <- bpvec(X=GRanges_annotSS_REF_strings,
                                          FUN=wmScore, object=wmAcceptorSites, BPPARAM=BPPARAM)

      # score ALT alleles for donor splice sites
      GRanges_annotSS_ALT_scores <- bpvec(X=GRanges_annotSS_ALT_strings,
                                          FUN=wmScore, object=wmAcceptorSites, BPPARAM=BPPARAM)

      # store the position of the allele respect to the position of the dinucleotide AG whose nucleotides occur at pos 1 and 2
      relposacceptor <- seq(-conservedPositions(wmAcceptorSites)[1]+1, -1, by=1)
      relposacceptor <- c(relposacceptor, seq(1, width(wmAcceptorSites)-length(relposacceptor)))

      SCOREss <- DataFrame(SCORE5ssREF=rep(NA_real_, length(GRanges_annotSS)),
                           SCORE5ssALT=rep(NA_real_, length(GRanges_annotSS)),
                           SCORE5ssPOS=rep(NA_integer_, length(GRanges_annotSS)),
                           SCORE3ssREF=round(GRanges_annotSS_REF_scores, digits=2),
                           SCORE3ssALT=round(GRanges_annotSS_ALT_scores, digits=2),
                           SCORE3ssPOS=relposacceptor[GRanges_annotSS$POS])

      ## incorporate the cryptic splice site annotations on synonymous variants into 'variantsVR'
      dummyDF[ssSNVmask, ] <- SCOREss
    } else
      warning(sprintf("Width of the 3' splice site region (%d) is not equal to the width of the weight matrix for acceptor sites (%d).",
                      wregion, width(wmAcceptorSites)))
  }

  ## coding synonymous variants

  message("Scoring potential cryptic splice sites in coding synonymous variants")

  if (any(variantsVR$CONSEQUENCE %in% "synonymous")) {
    ## %in% avoids NAs when comparing with them (THIS MASK IS ALSO USED BELOW !!)
    synonymousSNVmask <- variantsVR$CONSEQUENCE %in% "synonymous"

    GRanges_SY <- variantsVR[synonymousSNVmask]

    # retrieve regions around the allele potentially involving cryptic donor sites
    wregion <- width(wmDonorSites)*2-1
    GRanges_SY_window_donor <- resize(GRanges_SY, width=width(GRanges_SY)+wregion-1, fix="center") 
    GRanges_SY_donor_strings <- getSeq(bsgenome, GRanges_SY_window_donor)
    
    # replace the variant by the alternate allele. This requires creating a DNAStringSetList,
    # from the varAllele column (DNAStringSet), which contains the ALT allele but strand adjusted
    GRanges_SY_donor_ALT_strings <- replaceAt(GRanges_SY_donor_strings,
                                              IRanges(width(wmDonorSites), width(wmDonorSites)),
                                              DNAStringSetList(strsplit(as.character(GRanges_SY_window_donor$varAllele), split="", fixed=TRUE)))
  
    # retrieve regions around the allele potentially involving cryptic acceptor sites
    wregion <- width(wmAcceptorSites)*2-1
    GRanges_SY_window_acceptor <- resize(GRanges_SY, width=width(GRanges_SY)+wregion-1, fix="center") 
    GRanges_SY_acceptor_strings <- getSeq(bsgenome, GRanges_SY_window_acceptor)

    GRanges_SY_acceptor_ALT_strings <- replaceAt(GRanges_SY_acceptor_strings,
                                                 IRanges(width(wmAcceptorSites), width(wmAcceptorSites)),
                                                 DNAStringSetList(strsplit(as.character(GRanges_SY_window_acceptor$varAllele), split="", fixed=TRUE)))
  
    # score synonymous ALT alleles for donor splice sites
    GRanges_SY_donor_ALT_scores <- bpvec(X=GRanges_SY_donor_ALT_strings,
                                         FUN=wmScore, object=wmDonorSites, BPPARAM=BPPARAM)
    GRanges_SY_donor_ALT_scores <- matrix(GRanges_SY_donor_ALT_scores, ncol=width(wmDonorSites), byrow=TRUE)
    GRanges_SY_donor_ALT_scores <- t(apply(GRanges_SY_donor_ALT_scores, 1, function(x) {
                                        maxsco <- maxpos <- NA_real_
                                        if (any(!is.na(x))) {
                                          maxpos <- which.max(x)
                                          maxsco <- x[maxpos]
                                        }
                                        c(maxsco, maxpos)
                                      }))

    # score synonymous ALT alleles for acceptor splice sites
    GRanges_SY_acceptor_ALT_scores <- bpvec(X=GRanges_SY_acceptor_ALT_strings,
                                            FUN=wmScore, object=wmAcceptorSites, BPPARAM=BPPARAM)
    GRanges_SY_acceptor_ALT_scores <- matrix(GRanges_SY_acceptor_ALT_scores, ncol=width(wmAcceptorSites), byrow=TRUE)
    GRanges_SY_acceptor_ALT_scores <- t(apply(GRanges_SY_acceptor_ALT_scores, 1, function(x) {
                                        maxsco <- maxpos <- NA_real_
                                        if (any(!is.na(x))) {
                                          maxpos <- which.max(x)
                                          maxsco <- x[maxpos]
                                        }
                                        c(maxsco, maxpos)
                                      }))

    # score synonymous REF alleles at the site of the highest score with the ALT allele
    posHighestSitesALT <- GRanges_SY_donor_ALT_scores[!is.na(GRanges_SY_donor_ALT_scores[, 2]), 2]
    GRanges_SY_donor_strings_at_ALT <- subseq(GRanges_SY_donor_strings[!is.na(GRanges_SY_donor_ALT_scores[, 2])],
                                                start=posHighestSitesALT, width=width(wmDonorSites))
    GRanges_SY_donor_REF_scores <- bpvec(X=GRanges_SY_donor_strings_at_ALT,
                                         FUN=wmScore, object=wmDonorSites, BPPARAM=BPPARAM)
  
    posHighestSitesALT <- GRanges_SY_acceptor_ALT_scores[!is.na(GRanges_SY_acceptor_ALT_scores[, 2]), 2]
    GRanges_SY_acceptor_strings_at_ALT <- subseq(GRanges_SY_acceptor_strings[!is.na(GRanges_SY_acceptor_ALT_scores[, 2])],
                                                   start=posHighestSitesALT, width=width(wmAcceptorSites))
    GRanges_SY_acceptor_REF_scores <- bpvec(X=GRanges_SY_acceptor_strings_at_ALT,
                                            FUN=wmScore, object=wmAcceptorSites, BPPARAM=BPPARAM)

    # store the position of the allele respect to the position of the dinucleotide GT whose nucleotides occur at pos 1 and 2
    relposdonor <- seq(-conservedPositions(wmDonorSites)[1]+1, -1, by=1)
    relposdonor <- c(relposdonor, seq(1, width(wmDonorSites)-length(relposdonor)))
    # store the position of the allele respect to the position of the dinucleotide AG whose nucleotides occur at pos 1 and 2
    relposacceptor <- seq(-conservedPositions(wmAcceptorSites)[1]+1, -1, by=1)
    relposacceptor <- c(relposacceptor, seq(1, width(wmAcceptorSites)-length(relposacceptor)))

    SCOREss_syn <- DataFrame(SCORE5ssREF=rep(NA, nrow(GRanges_SY_donor_ALT_scores)),
                             SCORE5ssALT=round(GRanges_SY_donor_ALT_scores[, 1], digits=2),
                             SCORE5ssPOS=relposdonor[width(wmDonorSites)-GRanges_SY_donor_ALT_scores[, 2]+1],
                             SCORE3ssREF=rep(NA, nrow(GRanges_SY_acceptor_ALT_scores)),
                             SCORE3ssALT=round(GRanges_SY_acceptor_ALT_scores[, 1], digits=2),
                             SCORE3ssPOS=relposacceptor[width(wmAcceptorSites)-GRanges_SY_acceptor_ALT_scores[, 2]+1])
    SCOREss_syn$SCORE5ssREF[!is.na(SCOREss_syn$SCORE5ssALT)] <- round(GRanges_SY_donor_REF_scores, digits=2)
    SCOREss_syn$SCORE3ssREF[!is.na(SCOREss_syn$SCORE3ssALT)] <- round(GRanges_SY_acceptor_REF_scores, digits=2)

    ## incorporate the cryptic splice site annotations on synonymous variants into 'variantsVR'
    dummyDF[synonymousSNVmask, ] <- SCOREss_syn
  }

  ## intronic variants

  message("Annotating potential cryptic splice sites in intronic variants")

  intronicSNVmask <- variantsVR$TYPE == "SNV" & variantsVR$LOCATION == "intron"
  if (any(intronicSNVmask)) {
    GRanges_intron_SNV <- variantsVR[intronicSNVmask] ## THIS MASK IS ALSO USED BELOW !!!

    ## adjust alternate allele for strand since the adjusted varAllele only exists for coding variants
    ## and the column ALT is not adjusted - adapted from VariantAnnotation/R/methods-predictCoding.R
    nstrand <- as.vector(strand(GRanges_intron_SNV) == "-")
    altAlleleStrandAdjusted <- DNAStringSetList(strsplit(alt(GRanges_intron_SNV), split="", fixed=TRUE))
    if (any(nstrand))
      altAlleleStrandAdjusted[nstrand] <- relist(complement(unlist(altAlleleStrandAdjusted[nstrand])),
                                                 altAlleleStrandAdjusted[nstrand])

    ## retrieve regions around the allele potentially involving cryptic donor sites
    wregion <- width(wmDonorSites)*2-1
    GRanges_intron_SNV_window_donor <- resize(GRanges_intron_SNV, width=width(GRanges_intron_SNV)+wregion-1, fix="center")
    GRanges_intron_SNV_donor_strings <- getSeq(bsgenome, GRanges_intron_SNV_window_donor)
    GRanges_intron_SNV_donor_ALT_strings <- replaceAt(GRanges_intron_SNV_donor_strings,
                                                      IRanges(width(wmDonorSites), width(wmDonorSites)),
                                                      altAlleleStrandAdjusted)

    ## retrieve regions around the allele potentially involving cryptic acceptor sites
    wregion <- width(wmAcceptorSites)*2-1
    GRanges_intron_SNV_window_acceptor <- resize(GRanges_intron_SNV, width=width(GRanges_intron_SNV)+wregion-1, fix="center") 
    GRanges_intron_SNV_acceptor_strings <- getSeq(bsgenome, GRanges_intron_SNV_window_acceptor)
    GRanges_intron_SNV_acceptor_ALT_strings <- replaceAt(GRanges_intron_SNV_acceptor_strings,
                                                         IRanges(width(wmAcceptorSites), width(wmAcceptorSites)),
                                                         altAlleleStrandAdjusted)

    ## score intronic ALT alleles for donor splice sites
    GRanges_intron_SNV_donor_ALT_scores <- bpvec(X=GRanges_intron_SNV_donor_ALT_strings,
                                                 FUN=wmScore, object=wmDonorSites, BPPARAM=BPPARAM)
    GRanges_intron_SNV_donor_ALT_scores <- matrix(GRanges_intron_SNV_donor_ALT_scores, ncol=width(wmDonorSites), byrow=TRUE)
    GRanges_intron_SNV_donor_ALT_scores <- t(apply(GRanges_intron_SNV_donor_ALT_scores, 1, function(x) {
                                                     maxsco <- maxpos <- NA_real_
                                                     if (any(!is.na(x))) {
                                                       maxpos <- which.max(x)
                                                       maxsco <- x[maxpos]
                                                     }
                                                     c(maxsco, maxpos)
                                                   }))

    ## score intronic ALT alleles for acceptor splice sites
    GRanges_intron_SNV_acceptor_ALT_scores <- bpvec(X=GRanges_intron_SNV_acceptor_ALT_strings,
                                                    FUN=wmScore, object=wmAcceptorSites, BPPARAM=BPPARAM)
    GRanges_intron_SNV_acceptor_ALT_scores <- matrix(GRanges_intron_SNV_acceptor_ALT_scores, ncol=width(wmAcceptorSites), byrow=TRUE)
    GRanges_intron_SNV_acceptor_ALT_scores <- t(apply(GRanges_intron_SNV_acceptor_ALT_scores, 1, function(x) {
                                                     maxsco <- maxpos <- NA_real_
                                                     if (any(!is.na(x))) {
                                                       maxpos <- which.max(x)
                                                       maxsco <- x[maxpos]
                                                     }
                                                     c(maxsco, maxpos)
                                                   }))
 
    ## score intronic REF alleles at the site of the highest score with the ALT allele
    posHighestSitesALT <- GRanges_intron_SNV_donor_ALT_scores[!is.na(GRanges_intron_SNV_donor_ALT_scores[, 2]), 2]
    GRanges_intron_SNV_donor_strings_at_ALT <- subseq(GRanges_intron_SNV_donor_strings[!is.na(GRanges_intron_SNV_donor_ALT_scores[, 2])],
                                                        start=posHighestSitesALT, width=width(wmDonorSites))
    GRanges_intron_SNV_donor_REF_scores <- bpvec(X=GRanges_intron_SNV_donor_strings_at_ALT,
                                                 FUN=wmScore, object=wmDonorSites, BPPARAM=BPPARAM)

    posHighestSitesALT <- GRanges_intron_SNV_acceptor_ALT_scores[!is.na(GRanges_intron_SNV_acceptor_ALT_scores[, 2]), 2]
    GRanges_intron_SNV_acceptor_strings_at_ALT <- subseq(GRanges_intron_SNV_acceptor_strings[!is.na(GRanges_intron_SNV_acceptor_ALT_scores[, 2])],
                                                           start=posHighestSitesALT, width=width(wmAcceptorSites))
    GRanges_intron_SNV_acceptor_REF_scores <- bpvec(X=GRanges_intron_SNV_acceptor_strings_at_ALT,
                                                    FUN=wmScore, object=wmAcceptorSites, BPPARAM=BPPARAM)

    ## store the position of the allele respect to the position of the dinucleotide GT whose nucleotides occur at pos 1 and 2
    relposdonor <- seq(-conservedPositions(wmDonorSites)[1]+1, -1, by=1)
    relposdonor <- c(relposdonor, seq(1, width(wmDonorSites)-length(relposdonor)))
    relposacceptor <- seq(-conservedPositions(wmAcceptorSites)[1]+1, -1, by=1)
    relposacceptor <- c(relposacceptor, seq(1, width(wmAcceptorSites)-length(relposacceptor)))

    SCOREss_intron_SNV <- DataFrame(SCORE5ssREF=rep(NA, nrow(GRanges_intron_SNV_donor_ALT_scores)),
                                    SCORE5ssALT=round(GRanges_intron_SNV_donor_ALT_scores[, 1], digits=2),
                                    SCORE5ssPOS=relposdonor[width(wmDonorSites)-GRanges_intron_SNV_donor_ALT_scores[, 2]+1],
                                    SCORE3ssREF=rep(NA, nrow(GRanges_intron_SNV_acceptor_ALT_scores)),
                                    SCORE3ssALT=round(GRanges_intron_SNV_acceptor_ALT_scores[, 1], digits=2),
                                    SCORE3ssPOS=relposacceptor[width(wmAcceptorSites)-GRanges_intron_SNV_acceptor_ALT_scores[, 2]+1])
    SCOREss_intron_SNV$SCORE5ssREF[!is.na(SCOREss_intron_SNV$SCORE5ssALT)] <- round(GRanges_intron_SNV_donor_REF_scores, digits=2)
    SCOREss_intron_SNV$SCORE3ssREF[!is.na(SCOREss_intron_SNV$SCORE3ssALT)] <- round(GRanges_intron_SNV_acceptor_REF_scores, digits=2)
  
    ## incorporate the cryptic splice site annotations on synonymous variants into 'variantsVR'
    dummyDF[intronicSNVmask, ] <- SCOREss_intron_SNV
  }

  dummyDF
}
