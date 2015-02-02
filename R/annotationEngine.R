## main function carrying out the annotation of genetic variants taken as an
## input GRanges object. It uses different functions from this package and
## from the VariantAnnotation package

annotationEngine <- function(variantsGR, param, BPPARAM=bpparam()) {
  
  if (length(variantsGR) == 0) {
    variantsGR_annotated <- variantsGR
    mcols(variantsGR_annotated) <- .emptyAnnotations()
    return(variantsGR_annotated)
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

  ##############################
  ##                          ##
  ## CLEAN UP VARIANT INFO    ##
  ##                          ##
  ##############################

  ## clean the variant information landed on the 'names' of the input variantsGR 'GenomicRanges'
  ## object to leave either dbSNP ids given to the input VCF in this field or a maximum of 20
  ## characters

  ## vnames <- names(variantsGR)
  ## vnames2 <- vnames
  ## vnames2 <- rep(NA_character_, length(vnames))

  ## mt <- gregexpr("^[a-zA-Z0-9]+:[0-9]+_[ACGT/]+", vnames)
  ## mtstart <- unlist(mt, use.names=FALSE)
  ## mtlength <- sapply(mt, attr, "match.length")
  ## vnames2[mtstart != -1] <- substr(vnames[mtstart != -1], mtstart[mtstart != -1], mtlength[mtstart != -1])

  ## mt <- gregexpr("RS=[a-zA-Z]+[0-9]+", vnames)
  ## mtstart <- sapply(mt, "[", 1)
  ## mtlength <- sapply(mt, attr, "match.length")
  ## mtlength <- sapply(mtlength, "[", 1)
  ## vnames2[mtstart != -1] <- substr(vnames[mtstart != -1], mtstart[mtstart != -1]+3, mtlength[mtstart != -1])

  ## mt <- gregexpr("[a-zA-Z]+[0-9]+", vnames)
  ## mtstart <- sapply(mt, "[", 1)
  ## mtlength <- sapply(mt, attr, "match.length")
  ## mtlength <- sapply(mtlength, "[", 1)
  ## vnames2[mtstart != -1] <- substr(vnames[mtstart != -1], mtstart[mtstart != -1], mtlength[mtstart != -1])

  ## wh <- nchar(vnames2) > 20
  ## vnames2[wh] <- paste0(substr(vnames2[wh], 1, 20), "...")

  ## names(variantsGR) <- vnames2

  ############################
  ##                        ##
  ## INDEX ORIGINAL VARIANT ##
  ##                        ##
  ############################

  mcols(variantsGR) <- cbind(mcols(variantsGR), DataFrame(IDX=1:length(variantsGR)))

  ##############################
  ##                          ##
  ## ANNOTATE TYPE OF VARIANT ##
  ##                          ##
  ##############################
  
  message("Annotating variant type (SNV, InDel, MNV)")
  mcols(variantsGR) <- cbind(mcols(variantsGR), typeOfVariants(variantsGR))

  ##############################
  ##                          ##
  ##  SNP-CENTRIC ANNOTATIONS ##
  ##                          ##
  ##############################

  ## do this before variants get replicated because of different functional annotations
  message(sprintf("Annotating dbSNP identifiers with %s", snpdb@data_pkgname))
  mcols(variantsGR) <- cbind(mcols(variantsGR),
                             annotateVariants(snpdb, variantsGR, param, BPPARAM=BPPARAM))

  #######################
  ##                   ##
  ## ANNOTATE LOCATION ##
  ##                   ##
  #######################

  ## at the moment we are not interested in intergenic variants and we also leave promoter region
  ## boundaries at their default value. This could be parametrized if needed by the 'VariantFilteringParam' input object
  message("Annotating location with VariantAnnotation::locateVariants()")
  ## variantsGR_annotated <- locateVariants(variantsGR, txdb, AllVariants(intergenic=IntergenicVariants(0, 0)))
  ## seqinfo(variantsGR_annotated) <- seqinfo(variantsGR) ## locateVariants() drops the SeqInfo information

  ## put back GRanges names which in general correspond to the variant identifier from the VCF file
  ## names(variantsGR_annotated) <- names(variantsGR)[variantsGR_annotated$QUERYID]

  ## put back metadata columns IDX, the REF and ALT alleles, the variant TYPE, the FILTER flag information and
  ## the dbSNP annotated identifier
  ## mcols(variantsGR_annotated) <- cbind(mcols(variantsGR_annotated),
  ##                                mcols(variantsGR[variantsGR_annotated$QUERYID])[, c("IDX", "REF", "ALT", "TYPE", "FILTER", "dbSNP")])
  located_variantsGR <- locateVariants(variantsGR, txdb, AllVariants(intergenic=IntergenicVariants(0, 0)))
  variantsGR_annotated <- variantsGR[located_variantsGR$QUERYID]
  variantsGR_annotated$LOCATION <- located_variantsGR$LOCATION
  variantsGR_annotated$LOCSTART <- located_variantsGR$LOCSTART
  variantsGR_annotated$QUERYID <- located_variantsGR$QUERYID
  variantsGR_annotated$TXID <- located_variantsGR$TXID
  variantsGR_annotated$CDSID <- located_variantsGR$CDSID
  variantsGR_annotated$GENEID <- located_variantsGR$GENEID

  ## if the argument 'allTranscripts' is set to 'FALSE' then keep only once identical variants
  ## annotated with the same type of region from the same gene but in different tx
  if (!allTranscripts) {
    selmcols <- c("QUERYID", "LOCATION", "GENEID")
    dupsmask <- duplicated(mcols(variantsGR_annotated)[, selmcols])
    variantsGR_annotated <- variantsGR_annotated[!dupsmask]
  }

  ## remove unnecessary (by now) metadata columns (LOCEND, QUERYID, PRECEDEID, FOLLOWID)
  ## selmcols <- c("IDX", "LOCATION", "LOCSTART", "TXID", "CDSID", "GENEID", "REF", "ALT", "TYPE", "FILTER", "dbSNP")
  ## mcols(variantsGR_annotated) <- mcols(variantsGR_annotated)[, selmcols]
  rmcols <- match("QUERYID", colnames(mcols(variantsGR_annotated)))
  mcols(variantsGR_annotated) <- mcols(variantsGR_annotated)[, -rmcols]

  ## annotate cDNA position where applicable
  maskexonic <- variantsGR_annotated$LOCATION %in% c("coding", "fiveUTR", "threeUTR")
  cDNAloc <- .cDNAloc(variantsGR_annotated[maskexonic], txdb)
  variantsGR_annotated$cDNALOC <- rep(NA_integer_, length(variantsGR_annotated))
  variantsGR_annotated$cDNALOC[maskexonic] <- cDNAloc

  ##############################
  ##                          ##
  ## ANNOTATE CODING VARIANTS ##
  ##                          ##
  ##############################
  
  message("Annotating coding variants VariantAnnotation::predictCoding()")
  ## following the help page of 'predictCoding()', expand to have one row per alternate allele
  variantsGR_annotated_coding <- VRanges()
  if (any(variantsGR_annotated$LOCATION == "coding")) {
    variantsGR_annotated_coding <- variantsGR_annotated[variantsGR_annotated$LOCATION == "coding"]
    ## eltlen <- elementLengths(variantsGR_annotated_coding$ALT)
    ## variantsGR_annotated_coding_exp <- rep(variantsGR_annotated_coding, eltlen)

    ## harmonize sequence levels between variants, annotations and genome sequence
    origVarLevelsStyle <- seqlevelsStyle(variantsGR_annotated_coding)
    seqlevelsStyle(variantsGR_annotated_coding) <- seqlevelsStyle(bsgenome)
    origTxDbLevelsStyle <- seqlevelsStyle(txdb)
    seqlevelsStyle(txdb) <- seqlevelsStyle(bsgenome)

    commonChr <- intersect(seqlevels(variantsGR_annotated_coding), seqlevels(bsgenome))
    if (any(genome(variantsGR_annotated_coding)[commonChr] != genome(bsgenome)[commonChr])) {
      warning(sprintf("Assumming %s represent the same genome build.",
                      paste(c(unique(genome(variantsGR_annotated_coding)[commonChr]), unique(genome(bsgenome)[commonChr])),
                            collapse=" and ")))
      ## this generates conflicts with genome(txdb) in predictCoding() but it works
      ## if we do not do it, so we just issue the warning
      ## genome(variantsGR_annotated_coding_exp) <- genome(bsgenome)
    }

    ## remove columns that are again created by 'predictCoding()' to avoid this function prompting an error
    ## selmcols <- c("IDX", "LOCATION", "LOCSTART", "cDNALOC", "REF", "ALT", "TYPE", "FILTER", "dbSNP")
    ## mcols(variantsGR_annotated_coding_exp) <- mcols(variantsGR_annotated_coding_exp)[, selmcols]
    ## GRanges_coding_uq <- predictCoding(variantsGR_annotated_coding_exp, txdb,
    ##                                    seqSource=bsgenome, varAllele=unlist(variantsGR_annotated_coding$ALT))

    variantsGR_annotated_coding_exp <- as(variantsGR_annotated_coding, "GRanges")
    rmcols <- match(c("TXID", "CDSID", "GENEID"), colnames(mcols(variantsGR_annotated_coding_exp)))
    mcols(variantsGR_annotated_coding_exp) <- mcols(variantsGR_annotated_coding_exp)[, -rmcols]

    GRanges_coding_uq <- predictCoding(query=variantsGR_annotated_coding_exp,
                                       subject=txdb, seqSource=bsgenome,
                                       varAllele=DNAStringSet(alt(variantsGR_annotated_coding)))

    seqlevelsStyle(GRanges_coding_uq) <- origVarLevelsStyle
    seqlevelsStyle(txdb) <- origTxDbLevelsStyle
  
    ## if the argument 'allTranscripts' is set to 'FALSE' then keep only once identical variants
    ## annotated with the same coding consequence to the same gene but in different tx
    if (!allTranscripts)
      GRanges_coding_uq <- GRanges_coding_uq[!duplicated(mcols(GRanges_coding_uq)[, c("QUERYID", "GENEID", "CONSEQUENCE")])]

    variantsGR_annotated_coding <- variantsGR_annotated_coding[GRanges_coding_uq$QUERYID]

    ## some TXID, CDSID and GENEID annotations are catched by
    ## VariantAnnotation::locateVariants() and not by VariantAnnotation::predictCoding()
    ## and the other way araound. let's make the union of those and prompt and internal error
    ## when both annotations disagree.
    mask <- !is.na(variantsGR_annotated_coding$GENEID) & !is.na(GRanges_coding_uq$GENEID)
    if (any(variantsGR_annotated_coding$GENEID[mask] != GRanges_coding_uq$GENEID[mask]))
      stop("Internal error: enes annotated by VariantAnnotation::locateVariants() and VariantAnnotation::predictCoding() are different\n")

    mask <- is.na(variantsGR_annotated_coding$TXID)
    variantsGR_annotated_coding$TXID[mask] <- GRanges_coding_uq$TXID[mask]
    mask <- is.na(variantsGR_annotated_coding$CDSID)
    variantsGR_annotated_coding$CDSID[mask] <- GRanges_coding_uq$CDSID[mask]
    mask <- is.na(variantsGR_annotated_coding$GENEID)
    variantsGR_annotated_coding$GENEID[mask] <- GRanges_coding_uq$GENEID[mask]

    ## add coding annotations
    variantsGR_annotated_coding$CDSLOC <- GRanges_coding_uq$CDSLOC
    variantsGR_annotated_coding$PROTEINLOC <- GRanges_coding_uq$PROTEINLOC
    variantsGR_annotated_coding$REFCODON <- GRanges_coding_uq$REFCODON
    variantsGR_annotated_coding$VARCODON <- GRanges_coding_uq$VARCODON
    variantsGR_annotated_coding$REFAA <- GRanges_coding_uq$REFAA
    variantsGR_annotated_coding$VARAA <- GRanges_coding_uq$VARAA
    variantsGR_annotated_coding$varAllele <- GRanges_coding_uq$varAllele
    variantsGR_annotated_coding$CONSEQUENCE <- GRanges_coding_uq$CONSEQUENCE
  
    ## annotate codon usage difference in synonymous mutations

    ## if there are coding synonymous variants, add two metadatacolumns: CUREF and CUALT 
 
    ## GRanges_coding_uq_codonusage <- GRanges_coding_uq
    ## mcols(GRanges_coding_uq_codonusage) <- cbind(mcols(GRanges_coding_uq_codonusage),
    ##                                              DataFrame(CUREF=rep(NA_real_, length(GRanges_coding_uq_codonusage)),           
    ##                                                        CUALT=rep(NA_real_, length(GRanges_coding_uq_codonusage))))
    variantsGR_annotated_coding$CUREF <- NA_real_
    variantsGR_annotated_coding$CUALT <- NA_real_

    if (any(variantsGR_annotated_coding$CONSEQUENCE == "synonymous")){
    
      message("Annotating codon usage frequencies in coding synonymous variants")
    
      ## GRanges_coding_uq_codonusage <- GRanges_coding_uq
    
      ## create a mask for those variants that are codying and synonymous

      ## mask <- GRanges_coding_uq_codonusage$CONSEQUENCE %in% "synonymous"
      mask <- variantsGR_annotated_coding$CONSEQUENCE %in% "synonymous"
    
      ## extract the reference codons and alt codons from the DNStringSet
   
      ## ref_codons <- unname(as.character(GRanges_coding_uq_codonusage$REFCODON[mask]))
      ## alt_codons <- unname(as.character(GRanges_coding_uq_codonusage$VARCODON[mask]))
      ref_codons <- unname(as.character(variantsGR_annotated_coding$REFCODON[mask]))
      alt_codons <- unname(as.character(variantsGR_annotated_coding$VARCODON[mask]))

      ## obtain the vector with the frequencies for the reference and alternative codons   
      ref_codons_numeric <- codonusageTable[ref_codons]
      alt_codons_numeric <- codonusageTable[alt_codons]
    
      ## push the numerical vector to its corresponding variant to compare on CUREF and CUALT columns
      ## GRanges_coding_uq_codonusage$CUREF[mask] <- unname(ref_codons_numeric)
      ## GRanges_coding_uq_codonusage$CUALT[mask] <- unname(alt_codons_numeric)
      variantsGR_annotated_coding$CUREF[mask] <- unname(ref_codons_numeric)
      variantsGR_annotated_coding$CUALT[mask] <- unname(alt_codons_numeric)
        
    }
  
  }

  ## consolidate the annotations on coding variants, obatined with 'VariantAnnotation::predictCoding()',
  ##  with the rest of non-coding variants in a single VRanges object, replacing the one named
  ## 'variantsGR_annotated'
  variantsGR_annotated_noncoding <- variantsGR_annotated[variantsGR_annotated$LOCATION != "coding"]
  n.noncoding <- length(variantsGR_annotated_noncoding)
  dummyDF <- DataFrame(CDSLOC=IRanges(start=rep(-1, n.noncoding), end=rep(-1, n.noncoding)),
                       PROTEINLOC=IntegerList(as.list(rep(NA, n.noncoding))),
                       REFCODON=DNAStringSet(rep("", n.noncoding)),
                       VARCODON=DNAStringSet(rep("", n.noncoding)),
                       REFAA=AAStringSet(rep("", n.noncoding)),
                       VARAA=AAStringSet(rep("", n.noncoding)),
                       varAllele=DNAStringSet(rep("", n.noncoding)),
                       CONSEQUENCE=factor(rep(NA, n.noncoding), levels=c("nonsynonymous", "synonymous")),
                       CUREF=rep(NA_real_,n.noncoding),
                       CUALT=rep(NA_real_,n.noncoding))
  
  mcols(variantsGR_annotated_noncoding) <- cbind(mcols(variantsGR_annotated_noncoding), dummyDF)
  ## mcols(variantsGR_annotated_noncoding) <- mcols(variantsGR_annotated_noncoding)[, colnames(mcols(GRanges_coding_uq))]

  if (length(variantsGR_annotated_coding) > 0) {
    #mcols(GRanges_coding_uq) <- mcols(GRanges_coding_uq)[, colnames(mcols(variantsGR_annotated_noncoding))]
    ## mcols(GRanges_coding_uq_codonusage) <- mcols(GRanges_coding_uq_codonusage)[, colnames(mcols(variantsGR_annotated_noncoding))]
    #variantsGR_annotated <- c(GRanges_coding_uq, variantsGR_annotated_noncoding)
    ## variantsGR_annotated <- c(GRanges_coding_uq_codonusage, variantsGR_annotated_noncoding)
    variantsGR_annotated <- sort(c(variantsGR_annotated_coding, variantsGR_annotated_noncoding))
  } else
    variantsGR_annotated <- sort(variantsGR_annotated_noncoding)

  
  #############################################################
  ##                                                         ##
  ## ANNOTATE SPLICE SITES IN SYNONYMOUS & INTRONIC VARIANTS ##
  ##                                                         ##
  #############################################################
  
  ## add metadata columns in 'variantsGR_annotated' for cryptic ss annotations
  ## this should be optional once the shiny app is aware about present/absent annotations
  dummyDF <- DataFrame(CRYP5ssREF=rep(NA_real_, length(variantsGR_annotated)),
                       CRYP5ssALT=rep(NA_real_, length(variantsGR_annotated)),
                       CRYP5ssPOS=rep(NA_real_, length(variantsGR_annotated)),
                       CRYP3ssREF=rep(NA_real_, length(variantsGR_annotated)),
                       CRYP3ssALT=rep(NA_real_, length(variantsGR_annotated)),
                       CRYP3ssPOS=rep(NA_real_, length(variantsGR_annotated)))

  if (length(spliceSiteMatrices) == 2)
    dummyDF <- .scoreSpliceSiteVariants(variantsGR_annotated, spliceSiteMatrices, bsgenome, BPPARAM)
  mcols(variantsGR_annotated) <- cbind(mcols(variantsGR_annotated), dummyDF)

  ###############################
  ##                           ##
  ##  GENE-CENTRIC ANNOTATIONS ##
  ##                           ##
  ###############################

  mcols(variantsGR_annotated) <- cbind(mcols(variantsGR_annotated),
                                       annotateVariants(orgdb, variantsGR_annotated, param))

  #####################################
  ##                                 ##
  ##  TRANSCRIPT-CENTRIC ANNOTATIONS ##
  ##                                 ##
  #####################################

  mcols(variantsGR_annotated) <- cbind(mcols(variantsGR_annotated),
                                       annotateVariants(txdb, variantsGR_annotated, param))

  ## TODO: REPLACE THIS BY THE HGVS NOMENCLATURE !!!!!
  ########################
  ##                    ##
  ##        CDS         ##
  ##                    ##
  ########################

  mcols(variantsGR_annotated) <- cbind(mcols(variantsGR_annotated), DataFrame(CDS=rep(NA_character_, length(variantsGR_annotated))))

  ## for coding variants we use the 'varAllele' column which is already adjusted for strand
  refAllele <- as.character(adjustForStrandSense(variantsGR_annotated[variantsGR_annotated$LOCATION == "coding"],
                                                 ref(variantsGR_annotated)[variantsGR_annotated$LOCATION == "coding"]))
  locAllele <- as.character(start(variantsGR_annotated$CDSLOC[variantsGR_annotated$LOCATION == "coding"]))
  altAllele <- as.character(variantsGR_annotated$varAllele[variantsGR_annotated$LOCATION == "coding"])
  variantsGR_annotated$CDS[variantsGR_annotated$LOCATION == "coding"] <- paste0(refAllele, locAllele, altAllele)

  ## for non-coding variants we have to adjust for strand both, reference and alternative alleles
  ## THIS IS PROBABLY REDUNDANT AS THE VRanges CONSTRUCTOR ALREAD ADJUSTS FOR THIS (???)
  refAllele <- as.character(adjustForStrandSense(variantsGR_annotated[variantsGR_annotated$LOCATION != "coding"],
                                                 ref(variantsGR_annotated)[variantsGR_annotated$LOCATION != "coding"]))
  altAllele <- adjustForStrandSense(variantsGR_annotated[variantsGR_annotated$LOCATION != "coding"],
                                    alt(variantsGR_annotated)[variantsGR_annotated$LOCATION != "coding"])
  ## this is not necessary anymore as the VRanges is already expanded to variant x allele
  ## eltlen <- elementLengths(altAllele)
  ## altAlleleChr <- rep(NA_character_, length(altAllele))
  ## altAlleleChr[eltlen == 1] <- as.character(unlist(altAllele[eltlen == 1]))
  ## altAlleleChr[eltlen > 1] <- unlist(bplapply(altAllele[eltlen > 1], paste, collapse=", ", BPPARAM=BPPARAM),
  ##                                    use.names=FALSE)

  variantsGR_annotated$CDS[variantsGR_annotated$LOCATION != "coding"] <- paste(refAllele, "->", altAllele)

  ########################
  ##                    ##
  ##     AminoAcid      ##
  ##                    ##
  ########################

  mcols(variantsGR_annotated) <- cbind(mcols(variantsGR_annotated),
                                       aminoAcidChanges(variantsGR_annotated, radicalAAchangeMatrix))

  ########################
  ##                    ##
  ##  OTHER ANNOTATIONS ##
  ##                    ##
  ########################
  
  for (i in seq(along=otherAnnotations)) {
    message(sprintf("Annotating with %s", names(otherAnnotations)[i]))
    mcols(variantsGR_annotated) <- cbind(mcols(variantsGR_annotated),
                                         annotateVariants(otherAnnotations[[i]],
                                                          variantsGR_annotated,
                                                          param,
                                                          BPPARAM=BPPARAM))
  }

  ## this seems to be nicely transferred during the coercion to VRanges
  ## restore the seqinfo data
  ## seqinfo(variantsGR_annotated) <- seqinfo(variantsGR)

  return(variantsGR_annotated)
}



###########
## Annotate dbSNP identifiers
####

## an odd thing here is that getSNPlocs needs not the explicit SNPlocs object (annObj) as input argument
setMethod("annotateVariants", signature(annObj="SNPlocs"),
          function(annObj, variantsGR, param, BPPARAM=bpparam()) {
            if (!"TYPE" %in% colnames(mcols(variantsGR))) {
              stop("Variant type (SNV, InDel, MNV) has not been annotated.")
            }
            seqlevelsStyle(variantsGR) <- seqlevelsStyle(annObj)
            masksnp <- variantsGR$TYPE == "SNV"
            rsids_list <- loc2rsid(variantsGR[masksnp], BPPARAM=BPPARAM)
            rsids <- rep(NA_character_, times=length(variantsGR))
            elen <- elementLengths(rsids_list)
            rsids[masksnp][elen == 1] <- as.character(rsids_list[elen == 1])

            ## paste together multiple dbSNP identifiers
            rsids[masksnp][elen > 1] <- unlist(bplapply(rsids_list[elen > 1], paste, collapse=", ", BPPARAM=BPPARAM),
                                               use.names=FALSE)
            rsids[rsids == ""] <- NA_character_
            return(DataFrame(dbSNP=rsids))
          })

###########
## Annotate PolyPhen2 predictions
####

setMethod("annotateVariants", signature(annObj="PolyPhenDb"),
          function(annObj, variantsGR, param, coding=TRUE, BPPARAM=bpparam()) {
            PolyPhen2 <- rep(NA_character_, length(variantsGR))
            if (!coding)
              return(DataFrame(PolyPhen2=PolyPhen2))

            rsids <- unique(variantsGR$dbSNP)
            rsids <- rsids[!is.na(rsids)]
            if (length(rsids) > 0) {
              pp <- select(annObj, keys=rsids,
                           cols=c("TRAININGSET", "PREDICTION", "PPH2PROB"))
              mask <- pp$TRAININGSET == "humvar"
              pphumvar <- pp[mask, ]
              mt <- match(variantsGR$dbSNP, pphumvar$RSID)
              PolyPhen2[!is.na(mt)] <- pphumvar$PREDICTION[mt[!is.na(mt)]]
            }
            DataFrame(PolyPhen2=PolyPhen2)
          })

############
## Annotate PROVEAN predictions (former SIFT)
#####

setMethod("annotateVariants", signature(annObj="PROVEANDb"),
          function(annObj, variantsGR, param, coding=TRUE, BPPARAM=bpparam()) {
            PROVEAN <- rep(NA_character_, length(variantsGR))
            if (!coding)
              return(DataFrame(PROVEAN=PROVEAN))

            rsids <- unique(variantsGR$dbSNP)
            rsids <- rsids[!is.na(rsids)]
            if (length(rsids) > 0) {
              rsids <- gsub("rs", "", rsids)
              pv <- select(annObj, keys=rsids, columns="PROVEANPRED")
              mt <- match(gsub("rs", "", variantsGR$dbSNP), pv$DBSNPID)
              PROVEAN[!is.na(mt)] <- pv$PROVEANPRED[mt[!is.na(mt)]]
            }
            DataFrame(PROVEAN=PROVEAN)
          })

############
## Annotate MAF values
#####

## revise this according to the previous two methods!!!
setMethod("annotateVariants", signature(annObj="MafDb"),
          function(annObj, variantsGR, param, BPPARAM=bpparam()) {

            ## get the MAF columns
            mafCols <- knownVariantsMAFcols(annObj) ## assumes all MAF column names contain 'AF'

            mafValues <- matrix(NA, nrow=length(variantsGR), ncol=length(mafCols),
                                dimnames=list(NULL, mafCols))
            varIDs <- variantsGR$dbSNP     ## fetch first by the annotated dbSNP identifier
            uniqVarIDs <- unique(varIDs)
            uniqVarIDs <- uniqVarIDs[!is.na(uniqVarIDs)]
            if (length(varIDs) > 0) {
              uniqMAFvalues <- fetchKnownVariantsByID(annObj, uniqVarIDs)
              mt <- match(varIDs, uniqMAFvalues$varID)
              mafValues[!is.na(mt), ] <- as.matrix(uniqMAFvalues[mt[!is.na(mt)], mafCols, drop=FALSE])

              ## for missing entries then fetch by given identifier
              varIDs <- names(variantsGR)
              missingdbsnpIDs <- unique(c(varIDs[is.na(variantsGR$dbSNP)],
                                          varIDs[!is.na(match(variantsGR$dbSNP, uniqMAFvalues$varID[is.na(uniqMAFvalues$chrom)]))]))
              missingdbsnpIDs <- missingdbsnpIDs[!is.na(missingdbsnpIDs)]
              if (length(missingdbsnpIDs) > 0) {
                uniqMAFvalues <- fetchKnownVariantsByID(annObj, missingdbsnpIDs)
                mt <- match(varIDs, uniqMAFvalues$varID)
                mafValues[!is.na(mt), ] <- as.matrix(uniqMAFvalues[mt[!is.na(mt)], mafCols, drop=FALSE])
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
          function(annObj, variantsGR, param, BPPARAM=bpparam()) {

            genelevel_annot <- DataFrame(GENE=character(), OMIM=character())
            geneIDs <- variantsGR$GENEID
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
          function(annObj, variantsGR, param, BPPARAM=bpparam()) {

            txlevel_annot <- DataFrame(TXNAME=character())
            txIDs <- variantsGR$TXID
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
          function(annObj, variantsGR, param, BPPARAM=bpparam()) {

            sco <- scores(annObj, variantsGR)

            DataFrame(phastCons=sco)
          })

############
## Annotate gene phylostratum
#####

setMethod("annotateVariants", signature(annObj="GenePhylostrataDb"),
          function(annObj, variantsGR, param, BPPARAM=bpparam()) {

            gps <- genePhylostratum(annObj, variantsGR$GENEID)

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

typeOfVariants <- function(variantsGR) {

  type <- factor(levels=c("InDel", "MNV", "SNV"))
  if (length(variantsGR) > 0) {
    ## variants with multiple alternate alleles should have all
    ## alt alleles as SNV to be annotated as SNV
    ## wref <- nchar(variantsGR$REF)
    ## walt <- nchar(variantsGR$ALT)
    type <- factor(rep("SNV", times=length(variantsGR)), levels=c("InDel", "MNV", "SNV"))
    ## type[any(wref != walt)] <- "InDel"
    ## type[all(wref == walt) & wref != 1] <- "MNV"
    ## stopifnot(all(wref[type == "SNV"] == 1)) ## QC
    type[isIndel(variantsGR)] <- "InDel"
    type[isSubstitution(variantsGR) & !isSNV(variantsGR)] <- "MNV"
  }

  DataFrame(TYPE=type)
}

## adapted from http://permalink.gmane.org/gmane.science.biology.informatics.conductor/48456
loc2rsid <- function(locs, BPPARAM=bpparam()) {
  getSNPcountFun <- getSNPlocsFun <- NULL
  tryCatch({
    getSNPcountFun <- get("getSNPcount", mode="function")
  }, error=function(err) {
    stop("'getSNPcount() function not found. This is a problem related to loading the SNP-centric package.")
  })

  tryCatch({
    getSNPlocsFun <- get("getSNPlocs", mode="function")
  }, error=function(err) {
    stop("'getSNPlocs() function not found. This is a problem related to loading the SNP-centric package.")
  })

  if (!is(locs, "GRanges"))
    stop("'locs' must be a GRanges object")
  if (!all(width(locs) == 1L))
    stop("all ranges in 'locs' must be of width 1")
  common_seqlevels <- intersect(seqlevels(locs), names(getSNPcountFun()))
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
                               snplocs <- getSNPlocsFun(seqname, as.GRanges=TRUE)
                               hits <- findOverlaps(locs2, snplocs) ## findOverlaps on a GRanges faster than findMatches on a vector
                               ## hits <- findMatches(locs2, snplocs$loc)
                               if (length(hits) > 0) {
                                   rsids <- paste0("rs", snplocs$RefSNP_id[subjectHits(hits)])
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
aminoAcidChanges <- function(variantsGR, rAAch) {
  aachange <- aachangetype <- rep(NA_character_, length(variantsGR))

  ## annotate non-synonymous coding changes
  whnonsyn <- which(variantsGR$CONSEQUENCE == "nonsynonymous")
  locaa <- rep(NA_character_, length(length(whnonsyn)))
  refaa <- as.character(variantsGR[whnonsyn]$REFAA)
  altaa <- as.character(variantsGR[whnonsyn]$VARAA)
  elen <- elementLengths(variantsGR[whnonsyn]$PROTEINLOC)
  
  locaa[elen == 1] <- as.character(unlist(variantsGR[whnonsyn]$PROTEINLOC[elen == 1], use.names=FALSE))

  ## location of a multiple amino acid replacement is denoted by its position range
  locaa[elen > 1] <- sapply(variantsGR[whnonsyn]$PROTEINLOC[elen > 1],
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
  whsyn <- which(variantsGR$CONSEQUENCE == "synonymous")
  locaa <- rep(NA_character_, length(length(whsyn)))
  refaa <- as.character(variantsGR[whsyn]$REFAA)
  altaa <- as.character(variantsGR[whsyn]$VARAA)
  locaa <- as.character(unlist(variantsGR[whsyn]$PROTEINLOC, use.names=FALSE))
  aachange[whsyn] <- paste0(refaa, locaa, altaa)
  aachangetype[whsyn] <- "Conservative"

  ## annotate frameshift changes as "Radical"
  whframeshift <- which(variantsGR$CONSEQUENCE == "frameshift")
  locaa <- rep(NA_character_, length(length(whframeshift)))
  refaa <- as.character(variantsGR[whframeshift]$REFAA)
  altaa <- as.character(variantsGR[whframeshift]$VARAA)
  elen <- elementLengths(variantsGR[whframeshift]$PROTEINLOC)
  locaa[elen == 1] <- as.character(unlist(variantsGR[whframeshift]$PROTEINLOC[elen == 1], use.names=FALSE))
  ## location of a multiple amino acid replacement is denoted by its position range
  locaa[elen > 1] <- sapply(variantsGR[whframeshift]$PROTEINLOC[elen > 1],
                                      function(x) paste(range(x), collapse="-"))
  aachange[whframeshift] <- paste0(refaa, locaa, altaa)
  aachangetype[whframeshift] <- "Radical"

  ## annotate nonsense changes as "Radical"
  whnonsense <- which(variantsGR$CONSEQUENCE == "nonsense")
  locaa <- rep(NA_character_, length(length(whnonsense)))
  refaa <- as.character(variantsGR[whnonsense]$REFAA)
  altaa <- as.character(variantsGR[whnonsense]$VARAA)
  locaa <- as.character(unlist(variantsGR[whnonsense]$PROTEINLOC, use.names=FALSE))
  aachange[whnonsense] <- paste0(refaa, locaa, altaa)
  aachangetype[whnonsense] <- "Radical"

  DataFrame(AAchange=aachange, AAchangeType=aachangetype)
}

.emptyAnnotations <- function() {
  DataFrame(LOCATION=factor(),
            LOCSTART=integer(),
            TXID=integer(),
            CDSID=integer(),
            GENEID=character(),
            REF=DNAStringSet(),
            ALT=DNAStringSetList(),
            TYPE=character(),
            FILTER=character(),
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
            CRYP5ssREF=numeric(),
            CRYP5ssALT=numeric(),
            CRYP5ssPOS=numeric(),
            CRYP3ssREF=numeric(),
            CRYP3ssALT=numeric(),
            CRYP3ssPOS=numeric(),
            GENE=character(),
            OMIM=character(),
            TXNAME=character(),
            CDS=character(),
            AAchange=character(),
            AAchangeType=character())
}

.scoreSpliceSiteVariants <- function(variantsGR, spliceSiteMatrices, bsgenome, BPPARAM=bpparam()) {

  ## adapt to sequence style and genome version from the input
  ## BSgenome object, thus assuming positions are based on the same
  ## genome even though might be differently specified (i.e., hg38 vs GRCh38)
  seqlevelsStyle(variantsGR) <- seqlevelsStyle(bsgenome)
  commonChr <- intersect(seqlevels(variantsGR), seqlevels(bsgenome))
  if (any(is.na(genome(variantsGR)))) {
    warning(sprintf("Assuming the genome build of the input variants is %s.", unique(genome(bsgenome)[commonChr])))
    genome(variantsGR) <- genome(bsgenome)
  } else if (any(genome(variantsGR)[commonChr] != genome(bsgenome)[commonChr])) {
    warning(sprintf("Assumming %s represent the same genome build.",
                    paste(c(unique(genome(variantsGR)[commonChr]), unique(genome(bsgenome)[commonChr])),
                          collapse=" and ")))
    genome(variantsGR) <- genome(bsgenome)
  }

  ## add metadata columns in 'variantsGR' for cryptic ss annotations
  dummyDF <- DataFrame(CRYP5ssREF=rep(NA_real_, length(variantsGR)),
                       CRYP5ssALT=rep(NA_real_, length(variantsGR)),
                       CRYP5ssPOS=rep(NA_real_, length(variantsGR)),
                       CRYP3ssREF=rep(NA_real_, length(variantsGR)),
                       CRYP3ssALT=rep(NA_real_, length(variantsGR)),
                       CRYP3ssPOS=rep(NA_real_, length(variantsGR)))

  wmDonorSites <- spliceSiteMatrices$wmDonorSites
  wmAcceptorSites <- spliceSiteMatrices$wmAcceptorSites

  ## coding synonymous variants

  message("Annotating potential cryptic splice sites in coding synonymous variants")
  eltlen <- elementLengths(variantsGR$ALT) ## DISCARD MULTIPLE ALT ALLELES BY NOW (SHOULD DO SOMETHING DIFFERENT)

  if (any(variantsGR$CONSEQUENCE %in% "synonymous" & eltlen == 1)) {
    ## %in% avoids NAs when comparing with them (THIS MASK IS ALSO USED BELOW !!)
    synonymousSNVmask <- variantsGR$CONSEQUENCE %in% "synonymous" & eltlen == 1

    GRanges_SY <- variantsGR[synonymousSNVmask]

    # retrieve regions around the allele potentially involving cryptic donor sites
    wregion <- width(wmDonorSites)*2-1
    GRanges_SY_window_donor <- resize(GRanges_SY, width=width(GRanges_SY)+wregion-1, fix="center") 
    GRanges_SY_donor_strings <- getSeq(bsgenome, GRanges_SY_window_donor)
    
    # So here we do the same but creating a DNAStringSetList, from the varAllele column (DNAStringSet), which contains the ALT allele but strand adjusted
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
    relposacceptor <- seq(-conservedPositions(wmAcceptorSites)[1]+1, -1, by=1)
    relposacceptor <- c(relposacceptor, seq(1, width(wmAcceptorSites)-length(relposacceptor)))

    CRYPss_syn <- DataFrame(CRYP5ssREF=rep(NA, nrow(GRanges_SY_donor_ALT_scores)),
                            CRYP5ssALT=round(GRanges_SY_donor_ALT_scores[, 1], digits=2),
                            CRYP5ssPOS=relposdonor[width(wmDonorSites)-GRanges_SY_donor_ALT_scores[, 2]+1],
                            CRYP3ssREF=rep(NA, nrow(GRanges_SY_acceptor_ALT_scores)),
                            CRYP3ssALT=round(GRanges_SY_acceptor_ALT_scores[, 1], digits=2),
                            CRYP3ssPOS=relposacceptor[width(wmAcceptorSites)-GRanges_SY_acceptor_ALT_scores[, 2]+1])
    CRYPss_syn$CRYP5ssREF[!is.na(CRYPss_syn$CRYP5ssALT)] <- round(GRanges_SY_donor_REF_scores, digits=2)
    CRYPss_syn$CRYP3ssREF[!is.na(CRYPss_syn$CRYP3ssALT)] <- round(GRanges_SY_acceptor_REF_scores, digits=2)

    ## incorporate the cryptic splice site annotations on synonymous variants into 'variantsGR'
    dummyDF[synonymousSNVmask, ] <- CRYPss_syn
  }

  ## intronic variants

  message("Annotating potential cryptic splice sites in intronic variants")
  eltlen <- elementLengths(variantsGR$ALT) ## DISCARD MULTIPLE ALT ALLELES BY NOW (SHOULD DO SOMETHING DIFFERENT)

  intronicSNVmask <- variantsGR$TYPE == "SNV" & variantsGR$LOCATION == "intron" & eltlen == 1
  if (any(intronicSNVmask)) {
    GRanges_intron_SNV <- variantsGR[intronicSNVmask] ## THIS MASK IS ALSO USED BELOW !!!

    ## adjust alternate allele for strand since the adjusted varAllele only exists for coding variants
    ## and the column ALT is not adjusted - adapted from VariantAnnotation/R/methods-predictCoding.R
    nstrand <- as.vector(strand(GRanges_intron_SNV) == "-")
    altAlleleStrandAdjusted <- GRanges_intron_SNV$ALT
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

    CRYPss_intron_SNV <- DataFrame(CRYP5ssREF=rep(NA, nrow(GRanges_intron_SNV_donor_ALT_scores)),
                                   CRYP5ssALT=round(GRanges_intron_SNV_donor_ALT_scores[, 1], digits=2),
                                   CRYP5ssPOS=relposdonor[width(wmDonorSites)-GRanges_intron_SNV_donor_ALT_scores[, 2]+1],
                                   CRYP3ssREF=rep(NA, nrow(GRanges_intron_SNV_acceptor_ALT_scores)),
                                   CRYP3ssALT=round(GRanges_intron_SNV_acceptor_ALT_scores[, 1], digits=2),
                                   CRYP3ssPOS=relposacceptor[width(wmAcceptorSites)-GRanges_intron_SNV_acceptor_ALT_scores[, 2]+1])
    CRYPss_intron_SNV$CRYP5ssREF[!is.na(CRYPss_intron_SNV$CRYP5ssALT)] <- round(GRanges_intron_SNV_donor_REF_scores, digits=2)
    CRYPss_intron_SNV$CRYP3ssREF[!is.na(CRYPss_intron_SNV$CRYP3ssALT)] <- round(GRanges_intron_SNV_acceptor_REF_scores, digits=2)
  
    ## incorporate the cryptic splice site annotations on synonymous variants into 'variantsGR'
    dummyDF[intronicSNVmask, ] <- CRYPss_intron_SNV
  }

  dummyDF
}
