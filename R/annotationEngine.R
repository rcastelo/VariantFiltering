## main function carrying out the annotation of genetic variants taken as an
## input GRanges object. It uses different functions from this package and
## from the VariantAnnotation package

annotationEngine <- function(variantsVR, param, cache=new.env(parent=emptyenv()),
                             BPPARAM=bpparam("SerialParam")) {
  
  if (length(variantsVR) == 0) {
    variantsVR_annotated <- variantsVR
    mcols(variantsVR_annotated) <- .emptyAnnotations(param)
    return(variantsVR_annotated)
  }

  bsgenome <- get(param$bsgenome)
  orgdb <- get(param$orgdb)
  txdb <- get(param$txdb)
  snpdb <- sapply(param$snpdb, function(pkg) get(pkg))
  weightMatrices <- param$weightMatrices
  radicalAAchangeMatrix <- param$radicalAAchangeMatrix
  allTranscripts <- param$allTranscripts
  otherAnnotations <- sapply(param$otherAnnotations, function(pkg) get(pkg))
  codonusageTable <- param$codonusageTable
  geneticCode <- param$geneticCode

  ###########################
  ##                       ##
  ## CLEAN UP VARIANT INFO ##
  ##                       ##
  ###########################

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

    mt <- gregexpr("^[0-9]+$", vnames)
    mtstart <- sapply(mt, "[", 1)
    mtlength <- sapply(mt, attr, "match.length")
    mtlength <- sapply(mtlength, "[", 1)
    vnames2[mtstart != -1] <- substr(vnames[mtstart != -1], mtstart[mtstart != -1], mtlength[mtstart != -1])

    wh <- nchar(vnames2) > 20
    vnames2[wh] <- paste0(substr(vnames2[wh], 1, 20), "...")

    variantsVR$VARID <- vnames2
  }
  mcols(mcols(variantsVR)) <- DataFrame(TAB=rep(NA, ncol(mcols(variantsVR))))
  mt <- match("VARID", colnames(mcols(variantsVR)))
  stopifnot(all(!is.na(mt))) ## QC
  mcols(mcols(variantsVR))$TAB[mt] <- rep("Genome", length(mt))

  ## initialize annotation metadata
  annotationmetadata <- list(filters=list(),
                             cutoffs=CutoffsList(list()),
                             sortings=CutoffsList(criterion=factor("position", levels="position"),
                                                  decreasing=FALSE))

  ##############################
  ##                          ##
  ## ANNOTATE TYPE OF VARIANT ##
  ##                          ##
  ##############################
  
  message("Annotating variant type (SNV, Insertion, Deletion, MNV, Delins)")
  tovannot <- typeOfVariants(variantsVR)
  mcols(variantsVR) <- cbind(mcols(variantsVR), tovannot)
  if (!is.null(metadata(tovannot)$filters))
    annotationmetadata$filters <- c(annotationmetadata$filters, metadata(tovannot)$filters)
  if (!is.null(metadata(tovannot)$cutoffs))
    annotationmetadata$cutoffs <- c(annotationmetadata$cutoffs, metadata(tovannot)$cutoffs)

  #############################
  ##                         ##
  ## SNP-CENTRIC ANNOTATIONS ##
  ##                         ##
  #############################

  ## do this before variants get replicated because of different location annotations
  ## TODO: try to do all these below within annotateVariants()
  dbSNPannot <- DataFrame(dbSNP=rep(NA_character_, times=length(variantsVR)))
  for (i in seq_len(length(snpdb))) {
    message(sprintf("Annotating dbSNP identifiers with %s", names(snpdb)[i]))
    res <- annotateVariants(snpdb[[i]], variantsVR, param, BPPARAM=BPPARAM)
    maskNAdbsnp <- is.na(res$dbSNP)
    maskNAannotdbsnp <- is.na(dbSNPannot$dbSNP)
    dbSNPannot$dbSNP[maskNAannotdbsnp & !maskNAdbsnp] <- res$dbSNP[maskNAannotdbsnp & !maskNAdbsnp]
    dbSNPannot$dbSNP[!maskNAannotdbsnp & !maskNAdbsnp] <-
      paste(dbSNPannot$dbSNP[!maskNAannotdbsnp & !maskNAdbsnp], res$dbSNP[!maskNAannotdbsnp & !maskNAdbsnp], sep=", ")
  }
  mcols(dbSNPannot) <- DataFrame(TAB="Genome")
  mcols(variantsVR) <- cbind(mcols(variantsVR), dbSNPannot)
  rm(dbSNPannot)
  dbsnpfilter <- function(x) !is.na(VariantFiltering::allVariants(x, groupBy="nothing")$dbSNP)
  attr(dbsnpfilter, "description") <- "Presence in dbSNP"
  attr(dbsnpfilter, "TAB") <- "Genome"
  environment(dbsnpfilter) <- baseenv()
  annotationmetadata$filters <- c(annotationmetadata$filters, list(dbSNP=dbsnpfilter))

  #######################
  ##                   ##
  ## ANNOTATE LOCATION ##
  ##                   ##
  #######################

  ## at the moment we are not interested in intergenic variants and we also leave promoter region
  ## boundaries at their default value. This could be parametrized if needed by the 'VariantFilteringParam' input object
  message("Annotating location with VariantAnnotation::locateVariants()")
  located_variantsGR <- .locateAllVariants(vfParam=param, query=as(variantsVR, "GRanges"),
                                           subject=txdb, cache=cache, BPPARAM=BPPARAM)
  variantsVR_annotated <- variantsVR[located_variantsGR$QUERYID] ## REPLACE variantsVR_annotated by variantsVR ??
  dtf <- DataFrame(LOCATION=located_variantsGR$LOCATION,
                   LOCSTART=located_variantsGR$LOCSTART,
                   LOCEND=located_variantsGR$LOCEND,
                   LOCSTRAND=strand(located_variantsGR),
                   QUERYID=located_variantsGR$QUERYID,
                   TXID=located_variantsGR$TXID,
                   CDSID=located_variantsGR$CDSID,
                   GENEID=located_variantsGR$GENEID)
  rm(located_variantsGR)
  mcols(dtf) <- DataFrame(TAB=c("Transcript", "Transcript", "Transcript", "Transcript", NA, "Transcript", "Transcript", "Gene"))
  locfilter <- function(x) {
                 allowedlocations <- names(VariantFiltering::cutoffs(x)$location)[VariantFiltering::cutoffs(x)$location]
                 VariantFiltering::allVariants(x, groupBy="nothing")$LOCATION %in% allowedlocations
               }
  attr(locfilter, "description") <- "Variant location (intergenic, promoter, coding, intron, etc.)"
  attr(locfilter, "TAB") <- "Transcript"
  environment(locfilter) <- baseenv()
  locmask <- do.call("names<-", list(rep(TRUE, nlevels(dtf$LOCATION)), levels(dtf$LOCATION)))
  metadata(dtf) <- list(filters=list(location=locfilter),
                        cutoffs=CutoffsList(location=locmask))
  mcols(variantsVR_annotated) <- cbind(mcols(variantsVR_annotated), dtf)
  if (!is.null(metadata(dtf)$filters))
    annotationmetadata$filters <- c(annotationmetadata$filters, metadata(dtf)$filters)
  if (!is.null(metadata(dtf)$cutoffs))
    annotationmetadata$cutoffs <- c(annotationmetadata$cutoffs, metadata(dtf)$cutoffs)
  rm(dtf)

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
  dtf <- DataFrame(cDNALOC=rep(NA_integer_, length(variantsVR_annotated)))
  dtf$cDNALOC[maskexonic] <- cDNAloc
  mcols(dtf) <- DataFrame(TAB="Transcript")
  mcols(variantsVR_annotated) <- cbind(mcols(variantsVR_annotated), dtf)
  rm(dtf)

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

    ## it is necessary to explicitly strand the annotated variants when calling predictCoding()
    ## to handle properly SNPs overlapping transcripts in opposite sense (e.g., rs549311269)
    strand(variantsVR_annotated_coding_exp) <- variantsVR_annotated_coding_exp$LOCSTRAND

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
    variantsVR_annotated_coding$CDSSTART <- rep(NA_integer_, length(variantsVR_annotated_coding))
    variantsVR_annotated_coding$CDSEND <- rep(NA_integer_, length(variantsVR_annotated_coding))
    uniqTxIDs <- unique(as.character(unlist(variantsVR_annotated_coding$TXID, use.names=FALSE)))
    tryCatch({
      cdsinfo <- suppressMessages(select(txdb, keys=uniqTxIDs, columns=c("CDSSTART", "CDSEND"), keytype="TXID"))
      cdsinfoStart <- split(cdsinfo$CDSSTART, cdsinfo$TXID)
      cdsinfoStart <- sapply(cdsinfoStart, min)
      cdsinfoEnd <- split(cdsinfo$CDSEND, cdsinfo$TXID)
      cdsinfoEnd <- sapply(cdsinfoEnd, max)
      stopifnot(all(names(cdsinfoStart) == names(cdsinfoEnd))) ## QC
      cdsinfo <- data.frame(TXID=names(cdsinfoStart),
                            CDSSTART=as.integer(cdsinfoStart),
                            CDSEND=as.integer(cdsinfoEnd),
                            stringsAsFactors=FALSE)
      mt <- match(variantsVR_annotated_coding$TXID, cdsinfo$TXID)
      variantsVR_annotated_coding$CDSSTART <- cdsinfo$CDSSTART[mt]
      variantsVR_annotated_coding$CDSEND <- cdsinfo$CDSEND[mt]
    }, error=function(err) {
      misk <- ifelse(length(uniqTxIDs) > 3, sprintf("%s, ...", paste(head(uniqTxIDs, n=3), collapse=", ")),
                     paste(uniqTxIDs, collapse=", "))
      warning(sprintf("Could not find any TXIDs (%s) in the transcript-centric annotation package %s\n",
                      misk, txdb$packageName))
    })

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
                       CDSSTART=rep(NA_integer_, n.noncoding),
                       CDSEND=rep(NA_integer_, n.noncoding),
                       CUREF=rep(NA_real_,n.noncoding),
                       CUALT=rep(NA_real_,n.noncoding))
  
  mcols(variantsVR_annotated_noncoding) <- cbind(mcols(variantsVR_annotated_noncoding), dummyDF)

  if (length(variantsVR_annotated_coding) > 0) {
    variantsVR_annotated <- sort(c(variantsVR_annotated_coding, variantsVR_annotated_noncoding))
  } else
    variantsVR_annotated <- sort(variantsVR_annotated_noncoding)

  rm(variantsVR_annotated_coding)
  rm(variantsVR_annotated_noncoding)

  ## annotate start and end positions of TX (to determine pre-mRNA positions)
  ## FIXME: with intergenic variants
  variantsVR_annotated$TXSTART <- rep(NA, length(variantsVR_annotated))
  variantsVR_annotated$TXEND <- rep(NA, length(variantsVR_annotated))
  uniqTxIDs <- unique(as.character(unlist(variantsVR_annotated$TXID, use.names=FALSE)))
  tryCatch({
    txinfo <- suppressMessages(select(txdb, keys=uniqTxIDs, columns=c("TXSTART", "TXEND"), keytype="TXID"))
    mt <- match(as.character(unlist(variantsVR_annotated$TXID, use.names=FALSE)), txinfo$TXID)
    variantsVR_annotated$TXSTART <- txinfo$TXSTART[mt]
    variantsVR_annotated$TXEND <- txinfo$TXEND[mt]
  }, error=function(err) {
      misk <- ifelse(length(uniqTxIDs) > 3, sprintf("%s, ...", paste(head(uniqTxIDs, n=3), collapse=", ")),
                     paste(uniqTxIDs, collapse=", "))
      warning(sprintf("Could not find any TXIDs (%s) in the transcript-centric annotation package %s\n",
                      misk, txdb$packageName))
  })

  ## mt <- match(c("CDSLOC", "PROTEINLOC", "REFCODON", "VARCODON", "REFAA", "VARAA",
  ##               "varAllele", "CONSEQUENCE", "CDSSTART", "CDSEND", "CUREF", "CUALT"),
  ##             colnames(mcols(variantsVR_annotated)))
  mt <- match(c("CONSEQUENCE", "CUREF", "CUALT"), colnames(mcols(variantsVR_annotated)))
  stopifnot(all(!is.na(mt))) ## QC
  mcols(mcols(variantsVR_annotated))$TAB[mt] <- rep("Protein", length(mt))
  mt <- match(c("TXSTART", "TXEND"), colnames(mcols(variantsVR_annotated)))
  stopifnot(all(!is.na(mt))) ## QC
  mcols(mcols(variantsVR_annotated))$TAB[mt] <- rep("Transcript", length(mt))
  
  confilter <- function(x) {
                 allowedconseq <- names(VariantFiltering::cutoffs(x)$consequence)[VariantFiltering::cutoffs(x)$consequence]
                 VariantFiltering::allVariants(x, groupBy="nothing")$CONSEQUENCE %in% allowedconseq
               }
  attr(confilter, "description") <- "Variant consequence (frameshift, nonsense, nonsynonymous, synonymous)"
  attr(confilter, "TAB") <- "Protein"
  environment(confilter) <- baseenv()
  conmask <- do.call("names<-", list(rep(TRUE, nlevels(mcols(variantsVR_annotated)$CONSEQUENCE)),
                                         levels(mcols(variantsVR_annotated)$CONSEQUENCE)))
  conmetadata <- list(filters=list(consequence=confilter),
                      cutoffs=CutoffsList(consequence=conmask))
  annotationmetadata$filters <- c(annotationmetadata$filters, conmetadata$filters)
  annotationmetadata$cutoffs <- c(annotationmetadata$cutoffs, conmetadata$cutoffs)

  cufilter <- function(x) {
                cuFC <- VariantFiltering::allVariants(x, groupBy="nothing")$CUALT /
                        VariantFiltering::allVariants(x, groupBy="nothing")$CUREF
                mask <- abs(log2(cuFC)) >= log2(VariantFiltering::cutoffs(x)$codonusageFC[1])
                mask[is.na(mask)] <- FALSE
                mask
              }
  attr(cufilter, "description") <- "Minimum codon usage fold change between synonymous reference and alternative codons"
  attr(cufilter, "TAB") <- "Protein"
  environment(cufilter) <- baseenv()
  cumetadata <- list(filters=list(codonusageFC=cufilter),
                     cutoffs=CutoffsList(codonusageFC=1))
  annotationmetadata$filters <- c(annotationmetadata$filters, cumetadata$filters)
  annotationmetadata$cutoffs <- c(annotationmetadata$cutoffs, cumetadata$cutoffs)

  #############################################################
  ##                                                         ##
  ## ANNOTATE SPLICE SITES IN SYNONYMOUS & INTRONIC VARIANTS ##
  ##                                                         ##
  #############################################################
  
  if (length(weightMatrices) > 0) {
    dummyDF <- .scoreBindingSiteVariants(variantsVR_annotated, weightMatrices, bsgenome, BPPARAM)

    mcols(variantsVR_annotated) <- cbind(mcols(variantsVR_annotated), dummyDF)
  }

  ##############################
  ##                          ##
  ## GENE-CENTRIC ANNOTATIONS ##
  ##                          ##
  ##############################

  genecentricannot <- annotateVariants(orgdb, variantsVR_annotated, param)
  mcols(variantsVR_annotated) <- cbind(mcols(variantsVR_annotated), genecentricannot)
  if (!is.null(metadata(genecentricannot)$filters))
    annotationmetadata$filters <- c(annotationmetadata$filters, metadata(genecentricannot)$filters)
  if (!is.null(metadata(genecentricannot)$cutoffs))
    annotationmetadata$cutoffs <- c(annotationmetadata$cutoffs, metadata(genecentricannot)$cutoffs)
  rm(genecentricannot)

  ####################################
  ##                                ##
  ## TRANSCRIPT-CENTRIC ANNOTATIONS ##
  ##                                ##
  ####################################

  mcols(variantsVR_annotated) <- cbind(mcols(variantsVR_annotated),
                                       annotateVariants(txdb, variantsVR_annotated, param))

  #################################
  ##                             ##
  ## PROTEIN-CENTRIC ANNOTATIONS ##
  ##                             ##
  #################################

  aachangesannot <- aminoAcidChanges(variantsVR_annotated, radicalAAchangeMatrix)
  mcols(variantsVR_annotated) <- cbind(mcols(variantsVR_annotated), aachangesannot)
  if (!is.null(metadata(aachangesannot)$filters))
    annotationmetadata$filters <- c(annotationmetadata$filters, metadata(aachangesannot)$filters)
  if (!is.null(metadata(aachangesannot)$cutoffs))
    annotationmetadata$cutoffs <- c(annotationmetadata$cutoffs, metadata(aachangesannot)$cutoffs)
  rm(aachangesannot)

  ######################
  ##                  ##
  ## HGVS ANNOTATIONS ##
  ##                  ##
  ######################

  mcols(variantsVR_annotated) <- cbind(mcols(variantsVR_annotated), variantHGVS(variantsVR_annotated))
                                      
  #######################
  ##                   ##
  ## OTHER ANNOTATIONS ##
  ##                   ##
  #######################
  
  for (i in seq_along(otherAnnotations)) {
    message(sprintf("Annotating with %s", names(otherAnnotations)[i]))
    ann <- annotateVariants(otherAnnotations[[i]], variantsVR_annotated,
                            param, BPPARAM=BPPARAM)
    if (!is.null(metadata(ann)$filters))
      annotationmetadata$filters <- c(annotationmetadata$filters, metadata(ann)$filters)
    if (!is.null(metadata(ann)$cutoffs))
      annotationmetadata$cutoffs <- c(annotationmetadata$cutoffs, metadata(ann)$cutoffs)
    if (!is.null(metadata(ann)$sortings))
      annotationmetadata$sortings$criterion <- factor(annotationmetadata$sortings$criterion[1],
                                                     levels=c(as.character(annotationmetadata$sortings$criterion), metadata(ann)$sortings))
    mcols(variantsVR_annotated) <- cbind(mcols(variantsVR_annotated), ann)
  }

  metadata(mcols(variantsVR_annotated)) <- annotationmetadata

  return(variantsVR_annotated)
}



##
## Annotate dbSNP identifiers
##

setMethod("annotateVariants", signature(annObj="SNPlocs"),
          function(annObj, variantsVR, param, BPPARAM=bpparam("SerialParam")) {
            if (!"TYPE" %in% colnames(mcols(variantsVR))) {
              stop("Variant type (SNV, Insertion, Deletion, MNV, Delins) has not been annotated.")
            }
            ## adapt to sequence style and genome version from the input
            ## SNPlocs object, thus assuming positions are based on the same
            ## genome even though might be differently specified (i.e., hg19 vs GRCh37.p13 or hg38 vs GRCh38)
            seqlevelsStyle(variantsVR) <- seqlevelsStyle(annObj)
            commonChr <- intersect(seqlevels(variantsVR), seqlevels(annObj))
            if (any(is.na(genome(variantsVR)))) {
              warning(sprintf("Assuming the genome build of the input variants is the one of the SNPlocs package (%s).", unique(genome(annObj)[commonChr])))
            genome(variantsVR) <- genome(annObj)
            } else if (any(genome(variantsVR)[commonChr] != genome(annObj)[commonChr])) {
              warning(sprintf("Assumming %s represent the same genome build between variants and the SNPlocs package, respectively.",
                               paste(c(unique(genome(variantsVR)[commonChr]),
                                       unique(genome(annObj)[commonChr])),
                                     collapse=" and ")))
              genome(variantsVR) <- genome(annObj)
            }

            ## discard annotating variants on chromosomes with different lengths from the
            ## chromosomes in the annotation object to avoid errors with 'snpByOverlaps()'
            masksnp <- rep(TRUE, length(variantsVR))
            slenVR <- seqlengths(variantsVR)[commonChr]
            slenAnnObj <- seqlengths(annObj)[commonChr]
            if (any(slenVR != slenAnnObj)) {
              if (sum(slenVR != slenAnnObj) == 1) {
                warning(sprintf("Chromosome %s has different lengths between the input VCF and the SNPlocs package. Variants in this chromosome will not be annotated with the SNPlocs package", paste(commonChr[which(slenVR != slenAnnObj)], collapse=", ")))
              } else {
                warning(sprintf("Chromosomes %s have different lengths between the input VCF and the SNPlocs package. Variants in these chromosomes will not be annotated with the SNPlocs package", paste(commonChr[which(slenVR != slenAnnObj)], collapse=", ")))
              }
              if (sum(slenVR == slenAnnObj) == 0)
                stop("None of the chromosomes in the input VCF file has the same length as the chromosomes in the SNPlocs package. The genome reference sequence employed to generate the VCF file was probably different from the one in the SNPlocs package.")
              commonLenChr <- commonChr[slenVR == slenAnnObj]
              masksnp <- decode(seqnames(variantsVR) %in% commonLenChr)
            }

            masksnp <- masksnp & (variantsVR$TYPE == "SNV")
            rsids <- rep(NA_character_, times=length(variantsVR))
            if (any(masksnp)) {
              tmpVR <- variantsVR[masksnp]
              tmpVR <- keepSeqlevels(tmpVR, as.character(unique(seqnames(tmpVR)))) ## drop sequence levels not in use
              rsids_list <- .loc2dbSNPid(annObj, tmpVR, BPPARAM=BPPARAM)
              rm(tmpVR)
              elen <- elementNROWS(rsids_list)
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
            ## adapt to sequence style and genome version from the input
            ## SNPlocs object, thus assuming positions are based on the same
            ## genome even though might be differently specified (i.e., hg19 vs GRCh37.p13 or hg38 vs GRCh38)
            seqlevelsStyle(variantsVR) <- seqlevelsStyle(annObj)
            commonChr <- intersect(seqlevels(variantsVR), seqlevels(annObj))
            if (any(is.na(genome(variantsVR)))) {
              warning(sprintf("Assuming the genome build of the input variants is the one of the XtraSNPlocs package (%s).", unique(genome(annObj)[commonChr])))
            genome(variantsVR) <- genome(annObj)
            } else if (any(genome(variantsVR)[commonChr] != genome(annObj)[commonChr])) {
              warning(sprintf("Assumming %s represent the same genome build between variants and the XtraSNPlocs package, respectively.",
                               paste(c(unique(genome(variantsVR)[commonChr]),
                                       unique(genome(annObj)[commonChr])),
                                     collapse=" and ")))
              genome(variantsVR) <- genome(annObj)
            }

            ## discard annotating variants on chromosomes with different lengths from the
            ## chromosomes in the annotation object to avoid errors with 'snpByOverlaps()'
            maskxtrasnp <- rep(TRUE, length(variantsVR))
            slenVR <- seqlengths(variantsVR)[commonChr]
            slenAnnObj <- seqlengths(annObj)[commonChr]
            if (any(slenVR != slenAnnObj)) {
              if (sum(slenVR != slenAnnObj) == 1) {
                warning(sprintf("Chromosome %s has different lengths between the input VCF and the SNPlocs package. Variants in this chromosome will not be annotated with the XtraSNPlocs package", paste(commonChr[which(slenVR != slenAnnObj)], collapse=", ")))
              } else {
                warning(sprintf("Chromosomes %s have different lengths between the input VCF and the SNPlocs package. Variants in these chromosomes will not be annotated with the XtraSNPlocs package", paste(commonChr[which(slenVR != slenAnnObj)], collapse=", ")))
              }
              if (sum(slenVR == slenAnnObj) == 0)
                stop("None of the chromosomes in the input VCF file has the same length as the chromosomes in the XtraSNPlocs package. The genome reference sequence employed to generate the VCF file was probably different from the one in the XtraSNPlocs package.")
              commonLenChr <- commonChr[slenVR == slenAnnObj]
              maskxtrasnp <- decode(seqnames(variantsVR) %in% commonLenChr)
            }

            maskxtrasnp <- maskxtrasnp & (variantsVR$TYPE != "SNV")
            rsids <- rep(NA_character_, times=length(variantsVR))
            if (any(maskxtrasnp)) {
              tmpVR <- variantsVR[maskxtrasnp]
              tmpVR <- keepSeqlevels(tmpVR, as.character(unique(seqnames(tmpVR)))) ## drop sequence levels not in use
              rsids_list <- .loc2XtraSNPid(annObj, tmpVR, BPPARAM=BPPARAM)
              rm(tmpVR)
              elen <- elementNROWS(rsids_list)
              rsids[maskxtrasnp][elen == 1] <- as.character(rsids_list[elen == 1])

              ## paste together multiple dbSNP identifiers
              rsids[maskxtrasnp][elen > 1] <- unlist(bplapply(rsids_list[elen > 1], paste, collapse=", ", BPPARAM=BPPARAM),
                                                     use.names=FALSE)
              rsids[rsids == ""] <- NA_character_
            }

            dtf <- DataFrame(dbSNP=rsids)
            mcols(dtf) <- DataFrame(TAB="Genome")

            dtf
          })

##
## Annotate PolyPhen2 predictions
##

setMethod("annotateVariants", signature(annObj="PolyPhenDb"),
          function(annObj, variantsVR, param, coding=TRUE, BPPARAM=bpparam("SerialParam")) {
            PolyPhen2 <- rep(NA_character_, length(variantsVR))
            if (!coding)
              return(DataFrame(PolyPhen2=PolyPhen2))

            rsids <- unique(variantsVR$dbSNP)
            rsids <- rsids[!is.na(rsids)]
            if (length(rsids) > 0) {
              tryCatch({
                pp <- suppressMessages(select(annObj, keys=rsids,
                                       cols=c("TRAININGSET", "PREDICTION", "PPH2PROB")))
              }, error=function(err) {
                misk <- ifelse(length(rsids) > 3, sprintf("%s, ...", paste(head(rsids, n=3), collapse=", ")),
                               paste(rsids, collapse=", "))
                warning(sprintf("Could not find any rsIDs (%s) in the PolyPhenDb package.\n", misk))
              })
              mask <- pp$TRAININGSET == "humvar"
              pphumvar <- pp[mask, ]
              mt <- match(variantsVR$dbSNP, pphumvar$RSID)
              PolyPhen2[!is.na(mt)] <- pphumvar$PREDICTION[mt[!is.na(mt)]]
            }

            dtf <- DataFrame(PolyPhen2=PolyPhen2)
            mcols(dtf) <- DataFrame(TAB="Protein")

            dtf
          })

##
## Annotate PROVEAN predictions (former SIFT)
##

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
              tryCatch({
                pv <- suppressMessages(select(annObj, keys=rsids, columns="PROVEANPRED"))
              }, error=function(err) {
                misk <- ifelse(length(rsids) > 3, sprintf("%s, ...", paste(head(rsids, n=3), collapse=", ")),
                               paste(rsids, collapse=", "))
                warning(sprintf("Could not find any rsIDs (%s) in the PROVEANDb package.\n", misk))
              })
              mt <- match(gsub("rs", "", variantsVR$dbSNP), pv$DBSNPID)
              PROVEAN[!is.na(mt)] <- pv$PROVEANPRED[mt[!is.na(mt)]]
            }
            dtf <- DataFrame(PROVEAN=PROVEAN)
            mcols(dtf) <- DataFrame(TAB="Protein")

            dtf
          })

##
## Annotate organism-level gene-centric features
##

##
## Get HGNC gene symbol and other gene-centric annotations
##

setMethod("annotateVariants", signature(annObj="OrgDb"),
          function(annObj, variantsVR, param, BPPARAM=bpparam("SerialParam")) {

            defgenecols <- c("SYMBOL", "OMIM")
            defgenecols <- defgenecols[defgenecols %in% columns(annObj)]
            if (length(defgenecols) < 1 && !"SYMBOL" %in% defgenecols)
              stop("Currently the gene-centric annotation package should have at least a SYMBOL column.")

            genelevel_annot <- rep(list(rep(NA_character_, times=0)), length(defgenecols))
            names(genelevel_annot) <- defgenecols
            genelevel_annot <- do.call("DataFrame", genelevel_annot)
            geneIDs <- variantsVR$GENEID
            geneKeytype <- param$geneKeytype
            maskNAs <- is.na(geneIDs)
            if (length(geneIDs) > 0) {
              ## if input IDs are NAs output should also be NAs and avoid querying malformed keys
              genelevel_annot <- rep(list(rep(NA_character_, times=length(geneIDs))), length(defgenecols))
              names(genelevel_annot) <- defgenecols
              genelevel_annot <- do.call("DataFrame", genelevel_annot)
              if (sum(!maskNAs) > 0) {
                if (is.na(geneKeytype)) {
                  geneKeytype <- "ENTREZID"
                  if (substr(geneIDs[!maskNAs][1], 1, 3) == "ENS")
                    geneKeytype <- "ENSEMBL"
                  else if (substr(geneIDs[!maskNAs][1], 1, 3) %in% c("NM_", "NP_", "NR_", "XM_", "XP", "XR_", "YP_"))
                    geneKeytype <- "REFSEQ"
                }
                uniqIDs <- unique(geneIDs[!maskNAs])
                tryCatch({
                  res <- suppressMessages(select(annObj, keys=as.character(uniqIDs), columns=defgenecols, keytype=geneKeytype))
                  for (colname in defgenecols) {
                    colxgeneID <- sapply(split(res[[colname]], res[[geneKeytype]]),
                                         function(x) {
                                           if (all(!is.na(x))) x <- paste(unique(x), collapse=", ") ; x
                                         })
                    genelevel_annot[[colname]][!maskNAs] <- colxgeneID[geneIDs[!maskNAs]]
                  }
                }, error=function(err) {
                  misk <- ifelse(length(uniqIDs) > 3, sprintf("%s, ...", paste(head(uniqIDs, n=3), collapse=", ")),
                                 paste(uniqIDs, collapse=", "))
                  warning(sprintf("Could not retrieve any %s keys (%s) from the gene-centric annotation package %s\n",
                                  geneKeytype, misk, annObj$packageName))
                })

              }
            }
            genelevel_annot$GENE <- genelevel_annot$SYMBOL
            genelevel_annot$SYMBOL <- NULL
            mcols(genelevel_annot) <- DataFrame(TAB=rep("Gene", ncol(genelevel_annot)))
            if ("OMIM" %in% colnames(genelevel_annot)) {
              omimfilter <- function(x) !is.na(VariantFiltering::allVariants(x, groupBy="nothing")$OMIM)
              attr(omimfilter, "description") <- "Presence in OMIM"
              attr(omimfilter, "TAB") <- "Gene"
              environment(omimfilter) <- baseenv()
              metadata(genelevel_annot) <- list(filters=list(OMIM=omimfilter))
            }

            genelevel_annot
          })

setMethod("annotateVariants", signature(annObj="TxDb"),
          function(annObj, variantsVR, param, BPPARAM=bpparam("SerialParam")) {

            txlevel_annot <- DataFrame(TXNAME=character())
            txIDs <- as.character(variantsVR$TXID)
            maskNAs <- is.na(txIDs)
            if (length(txIDs) > 0) {
              ## if input IDs are NAs output should also be NAs and avoid querying malformed keys
              txlevel_annot <- DataFrame(TXNAME=rep(NA_character_, times=length(txIDs)))
              if (sum(!maskNAs) > 0) {
                uniqTxIDs <- unique(txIDs[!maskNAs])
                tryCatch({
                  res <- suppressMessages(select(annObj, keys=uniqTxIDs, columns="TXNAME", keytype="TXID"))
                  txnamextxID <- sapply(split(res$TXNAME, res$TXID),
                                         function(x) paste(unique(x), collapse=", "))
                  txlevel_annot[!maskNAs, ] <- DataFrame(TXNAME=txnamextxID[txIDs[!maskNAs]])
                }, error=function(err) {
                  misk <- ifelse(length(uniqTxIDs) > 3, sprintf("%s, ...", paste(head(uniqTxIDs, n=3), collapse=", ")),
                                 paste(uniqTxIDs, collapse=", "))
                  warning(sprintf("Could not retrieve any keys (%s) from the transcript-centric annotation package %s\n",
                                  misk, annObj$packageName))
                })
              }
            }
            mcols(txlevel_annot) <- DataFrame(TAB=rep("Transcript", ncol(txlevel_annot)))

            txlevel_annot
          })

##
## Annotate with GScores
##

setMethod("annotateVariants", signature(annObj="GScores"),
          function(annObj, variantsVR, param, BPPARAM=bpparam("SerialParam")) {

            if ("default" %in% populations(annObj))
              defaultPopulation(annObj) <- "default"

            pop <- defaultPopulation(annObj)
            if (pop != "default")
              pop <- populations(annObj)

            cnames <- gscoresTag(annObj)
            if (any(pop != "default"))
              cnames <- paste0(cnames, pop)

            scodtf <- DataFrame(as.data.frame(matrix(NA_real_,
                                                     nrow=length(variantsVR),
                                                     ncol=length(pop),
                                                     dimnames=list(NULL, cnames))))
            snvmask <- rep(TRUE, length(variantsVR))
            if (gscoresNonSNRs(annObj))
              snvmask <- isSNV(variantsVR)

            ## need to check for NAs in 'snvmask' due to '*' ALT alleles
            ## see https://samtools.github.io/hts-specs/VCFv4.2.pdf (1.4.1.5)
            ## and https://software.broadinstitute.org/gatk/documentation/article.php?id=6926
            ## and https://support.bioconductor.org/p/80212
            if (any(!is.na(snvmask) & snvmask))
              scodtf[!is.na(snvmask) & snvmask, ] <- score(annObj,
                                                           ranges=variantsVR[!is.na(snvmask) & snvmask],
                                                           pop=pop, type="snrs")

            if (any(!is.na(snvmask) & !snvmask))
              scodtf[!is.na(snvmask) & !snvmask, ] <- score(annObj,
                                                            ranges=variantsVR[!is.na(snvmask) & !snvmask],
                                                            pop=pop, type="nonsnrs")

            mcols(scodtf) <- DataFrame(TAB=gscoresGroup(annObj))

            fstr <- sprintf("gscfilter <- function(x) { opf <- get(as.character(VariantFiltering::cutoffs(x)$%s$op[1])) ; opf(VariantFiltering::allVariants(x, groupBy=\"nothing\")$%s, VariantFiltering::cutoffs(x)$%s$value[1]) }", gscoresTag(annObj), gscoresTag(annObj), gscoresTag(annObj))
            if (length(pop) > 1)
              fstr <- sprintf("gscfilter <- function(x) { cutoffval <- rep(NA_real_, length(x)) ; mask <- rep(TRUE, length(x)) ; popmask <- VariantFiltering::cutoffs(x)$%s$popmask ; if (any(popmask)) { cutoffval <- do.call(get(as.character(VariantFiltering::cutoffs(x)$%s$pfun)), c(as.list(S4Vectors::mcols(VariantFiltering::allVariants(x, groupBy=\"nothing\"))[, names(popmask[popmask]), drop=FALSE]), na.rm=TRUE)) ; cutoffval[is.na(cutoffval)] <- -Inf ; opf <- get(as.character(VariantFiltering::cutoffs(x)$%s$op[1])) ; mask <- opf(cutoffval, VariantFiltering::cutoffs(x)$%s$value) ; mask[is.na(mask)] <- FALSE } ; mask }", gscoresTag(annObj), gscoresTag(annObj), gscoresTag(annObj), gscoresTag(annObj))
            eval(parse(text=fstr))
            attr(gscfilter, "description") <- sprintf("%s genomic scores", type(annObj))
            attr(gscfilter, "TAB") <- gscoresGroup(annObj)
            environment(gscfilter) <- baseenv()
            eval(parse(text=sprintf("flt <- list(%s=gscfilter)", gscoresTag(annObj))))
            eval(parse(text=sprintf("ctf <- CutoffsList(%s=CutoffsList(op=factor(\">=\", levels=c(\"==\", \"!=\", \">\", \"<\", \">=\", \"<=\")), value=0))", gscoresTag(annObj))))
            if (length(pop) > 1) {
              eval(parse(text=sprintf("ctf <- CutoffsList(%s=CutoffsList(popmask=setNames(rep(TRUE, ncol(scodtf)), colnames(scodtf)), pfun=factor(\"pmax\", levels=c(\"pmax\", \"pmin\")), op=factor(\"<=\", levels=c(\"==\", \"!=\", \">\", \"<\", \">=\", \"<=\")), value=%s))", gscoresTag(annObj), "0.5"))) ## REPLACE 0.5!!!
            }
            metadata(scodtf) <- list(filters=get("flt"), cutoffs=get("ctf"), sortings=type(annObj))

            scodtf
          })

##
## Annotate MAF values
##

setMethod("annotateVariants", signature(annObj="MafDb"),
          function(annObj, variantsVR, param, BPPARAM=bpparam("SerialParam")) {

            ## adapt to sequence style and genome version from the input
            ## MafDb object, thus assuming positions are based on the same
            ## genome even though might be differently specified (i.e., hg19 vs GRCh37.p13 or hg38 vs GRCh38)
            seqlevelsStyle(variantsVR) <- seqlevelsStyle(annObj)
            commonChr <- intersect(seqlevels(variantsVR), seqlevels(annObj))
            if (any(is.na(genome(variantsVR)))) {
              warning(sprintf("Assuming the genome build of the input variants is the one of the MafDb package (%s).", unique(genome(annObj)[commonChr])))
            genome(variantsVR) <- genome(annObj)
            } else if (any(genome(variantsVR)[commonChr] != genome(annObj)[commonChr])) {
              warning(sprintf("Assumming %s represent the same genome build between variants and the MafDb package, respectively.",
                               paste(c(unique(genome(variantsVR)[commonChr]),
                                       unique(genome(annObj)[commonChr])),
                                     collapse=" and ")))
              genome(variantsVR) <- genome(annObj)
            }

            ## get the MAF columns
            mafCols <- populations(annObj)

            mafValues <- DataFrame(as.data.frame(matrix(NA_real_, nrow=length(variantsVR),
                                                        ncol=length(mafCols),
                                                        dimnames=list(NULL, mafCols))))
            snvmask <- isSNV(variantsVR)
            if (any(snvmask))
              mafValues[snvmask, ] <- mcols(gscores(annObj, variantsVR[snvmask], mafCols, type="snrs"))

            if (any(!snvmask))
              mafValues[!snvmask, ] <- mcols(gscores(annObj, variantsVR[!snvmask], mafCols, type="nonsnrs"))

            colnames(mafValues) <- paste0(colnames(mafValues), annObj$tag) ## tag MAF columns with their data source
            mcols(mafValues) <- DataFrame(TAB=rep("MAF", ncol(mafValues)))
            maffilter <- function(x) {
                           maxMAFannot <- rep(NA_real_, length(x))
                           mask <- rep(TRUE, length(x))
                           mafpopmask <- VariantFiltering::cutoffs(x)$maxMAF$popmask
                           if (any(mafpopmask)) {
                             mtNoMAF <- NULL
                             maxMAFannot <- do.call(pmax,
                                                    c(as.list(S4Vectors::mcols(VariantFiltering::allVariants(x, groupBy="nothing"))[, names(mafpopmask[mafpopmask]), drop=FALSE]), na.rm=TRUE))
                             maxMAFannot[is.na(maxMAFannot)] <- -Inf
                             mask <- maxMAFannot <= VariantFiltering::cutoffs(x)$maxMAF$maxvalue
                             mask[is.na(mask)] <- FALSE
                           }
                           mask
                         }
            attr(maffilter, "description") <- "Maximum minor allele frequency"
            attr(maffilter, "TAB") <- "MAF"
            environment(maffilter) <- baseenv()
            cnAF <- rep(TRUE, ncol(mafValues))
            names(cnAF) <- colnames(mafValues)
            metadata(mafValues) <- list(filters=list(maxMAF=maffilter),
                                        cutoffs=CutoffsList(maxMAF=CutoffsList(popmask=cnAF, maxvalue=0.5)))

            mafValues
          })

##
## Annotate gene phylostratum
##

setMethod("annotateVariants", signature(annObj="GenePhylostrataDb"),
          function(annObj, variantsVR, param, BPPARAM=bpparam("SerialParam")) {

            gps <- genePhylostratum(annObj, variantsVR$GENEID)

            dtf <- DataFrame(GenePhylostratumTaxID=gps$TaxID,
                             GenePhylostratumIndex=gps$OldestPhylostratum,
                             GenePhylostratum=gps$Description)
            mcols(dtf) <- DataFrame(TAB=c(NA, NA, "Gene"))
            gpsfilter <- function(x) {
              VariantFiltering::allVariants(x, groupBy="nothing")$GenePhylostratumIndex <= as.integer(VariantFiltering::cutoffs(x)$genePhyloStratum[1])
            }
            attr(gpsfilter, "description") <- "Gene phylostratum"
            attr(gpsfilter, "TAB") <- "Gene"
            environment(gpsfilter) <- baseenv()
            allowedgps <- VariantFiltering::genePhylostrata(annObj)$Description
            ctf <- CutoffsList(genePhyloStratum=factor(allowedgps[length(allowedgps)],
                                                       levels=allowedgps))
            metadata(dtf) <- list(filters=list(genePhyloStratum=gpsfilter),
                                  cutoffs=ctf)
            dtf
          })


##
## Read matrix of radical versus conservative amino acid substitutions
##

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
    type[as.vector(isInsertion(variantsVR))] <- "Insertion"
    type[as.vector(isDeletion(variantsVR))] <- "Deletion"
    type[as.vector(isSubstitution(variantsVR) & !isSNV(variantsVR))] <- "MNV"
    type[as.vector(isDelins(variantsVR))] <- "Delins"
  }

  dtf <- DataFrame(TYPE=type)
  mcols(dtf) <- DataFrame(TAB="Genome")
  tovfilter <- function(x) {
                 allowedtypes <- names(VariantFiltering::cutoffs(x)$variantType)[VariantFiltering::cutoffs(x)$variantType]
                 VariantFiltering::allVariants(x, groupBy="nothing")$TYPE %in% allowedtypes
               }
  attr(tovfilter, "description") <- "Type of variant (SVN, Insertion, Deletion, MNV, Delins)"
  attr(tovfilter, "TAB") <- "Genome"
  environment(tovfilter) <- baseenv()
  vtmask <- do.call("names<-", list(rep(TRUE, nlevels(type)), levels(type)))
  metadata(dtf) <- list(filters=list(variantType=tovfilter),
                        cutoffs=CutoffsList(variantType=vtmask))
  dtf
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

  refAllele <- .adjustForStrandSense(variantsVR, ref(variantsVR))
  altAllele <- as.character(variantsVR$varAllele)

  maskCoding <- variantsVR$LOCATION == "coding"
  altAllele[!maskCoding] <- as.character(.adjustForStrandSense(variantsVR[!maskCoding], alt(variantsVR)[!maskCoding]))

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
    stopifnot(!is.null(variantsVR$LOCSTRAND)) ## QC 
    locStartAllele[mask] <- ifelse(variantsVR$LOCSTRAND[mask] == "+",
                                   as.integer(start(variantsVR)[mask]) - variantsVR$TXSTART[mask] + 1L,
                                   variantsVR$TXEND[mask] - as.integer(end(variantsVR)[mask]) + 1L)
    locEndAllele[mask] <- ifelse(variantsVR$LOCSTRAND[mask] == "+",
                                 as.integer(end(variantsVR)[mask]) - variantsVR$TXSTART[mask] + 1L,
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
  
  dtf <- DataFrame(HGVSg=gDesc, HGVSc=cDesc, HGVSp=pDesc)
  mcols(dtf) <- DataFrame(TAB=c("Genome", "Transcript", "Protein"))
  dtf
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

.loc2dbSNPid <- function(obj, locs, BPPARAM=bpparam("SerialParam")) {

  if (!is(locs, "GRanges"))
    stop("'locs' must be a GRanges object")

  mcols(locs) <- NULL
  locs <- as(locs, "GRanges") ## when 'locs' is a 'VRanges' object

  common_seqlevels <- intersect(seqlevels(locs), names(snpcount(obj)))
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
                               
                               locs3 <- snpsByOverlaps(obj, locs2, minoverlap=1L)
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
  elen <- elementNROWS(variantsVR[whnonsyn]$PROTEINLOC)
  
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
  elen <- elementNROWS(variantsVR[whframeshift]$PROTEINLOC)
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

  dtf <- DataFrame(AAchange=aachange, AAchangeType=aachangetype)
  mcols(dtf) <- DataFrame(TAB=c("Protein", "Protein"))
  aafilter <- function(x) {
                mask <- rep(TRUE, length(x))
                if (!is.null(VariantFiltering::cutoffs(x)$aaChangeType)) {
                  aachangetype <- VariantFiltering::allVariants(x, groupBy="nothing")$AAchangeType
                  mask <- is.na(aachangetype) | (aachangetype %in% VariantFiltering::cutoffs(x)$aaChangeType[1])
                }
                mask
               }
  attr(aafilter, "description") <- "Type of amino acid change (Conservative, Radical)"
  attr(aafilter, "TAB") <- "Protein"
  environment(aafilter) <- baseenv()
  metadata(dtf) <- list(filters=list(aaChangeType=aafilter),
                        cutoffs=CutoffsList(aaChangeType=factor("Conservative",
                                                                levels=c("Conservative", "Radical"))))
  dtf
}

.emptyAnnotations <- function(vfParObj) {
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
            GENE=character(),
            OMIM=character(),
            TXNAME=character(),
            HGVSg=character(),
            HGVSc=character(),
            HGVSp=character(),
            AAchange=character(),
            AAchangeType=character())
}
