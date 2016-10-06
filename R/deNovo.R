setMethod("deNovo", signature(param="VariantFilteringParam"),
          function(param, svparam=ScanVcfParam(),
                   use=c("everything", "complete.obs", "all.obs"),
                   BPPARAM=bpparam("SerialParam")) {

  ## store call for reproducing it later
  callobj <- match.call()
  callstr <- gsub(".local", "deNovo", deparse(callobj))

  ## fetch necessary parameters
  vcfFiles <- param$vcfFiles
  ped <- param$pedFilename
  seqInfos <- param$seqInfos
  txdb <- param$txdb
  bsgenome <- param$bsgenome
  sampleNames <- param$sampleNames
 
  if (!exists(as.character(substitute(BPPARAM))))
    stop(sprintf("Parallel back-end function %s given in argument 'BPPARAM' does not exist in the current workspace. Either you did not write correctly the function name or you did not load the package 'BiocParallel'.", as.character(substitute(BPPARAM))))
  
  if (length(vcfFiles) > 1)
    stop("More than one input VCF file is currently not supported. Please either merge the VCF files into a single one with software such as vcftools or GATK, or do the variant calling simultaneously on all samples, or proceed analyzing each file separately.")
  else if (length(vcfFiles) < 1)
    stop("A minimum of 1 vcf file has to be provided")

  if (is.na(ped))
    stop("Please specify a PED file name when building the parameter object.")

  pedDf <- .readPEDfile(ped)

  if (sum(pedDf$Phenotype == 1) != 2 || sort(pedDf$Gender[pedDf$Phenotype == 1]) != 1:2)
    stop("Current 'de novo' analysis requires two unaffected parents and one or more affected children.")

  unaff <- pedDf[pedDf$Phenotype == 1, ]
  unaff_dad <- unaff[unaff$Gender == 1, "IndividualID"]
  unaff_mom <- unaff[unaff$Gender == 2, "IndividualID"]
  aff <- pedDf[pedDf$Phenotype == 2, ]
  
  if (aff$FatherID != unaff_dad || aff$MotherID != unaff_mom)
    stop("Current 'de novo' analysis requires two unaffected parents and one or more affected children.")

  annotationCache <- new.env() ## cache annotations when using VariantAnnotation::locateVariants()
  annotated_variants <- VRanges()
  open(vcfFiles[[1]])
  n.var <- 0
  while (nrow(vcf <- readVcf(vcfFiles[[1]], genome=seqInfos[[1]], param=svparam))) {

    ## insert an index for each variant in the VCF file
    info(header(vcf)) <- rbind(info(header(vcf)),
                               DataFrame(Number=1, Type="Integer",
                                         Description="Variant index in the VCF file.",
                                         row.names="VCFIDX"))
    info(vcf)$VCFIDX <- (n.var+1):(n.var+nrow(vcf))
    varIDs <- rownames(vcf)

    n.var <- n.var + nrow(vcf)
  
    mask <- .deNovoMask(vcf, pedDf, use)

    if (any(mask)) {

      ## filter out variants that do not segregate as a "de novo" trait
      vcf <- vcf[mask, ]
     
      ## coerce the VCF object to a VRanges object
      variants <- as(vcf, "VRanges")

      ## since the conversion of VCF to VRanges strips the VCF ID field, let's put it back
      variants$VARID <- varIDs[variants$VCFIDX]

      ## harmonize Seqinfo data between variants, annotations and reference genome
      variants <- .matchSeqinfo(variants, txdb, bsgenome)

      ## annotate variants
      annotated_variants <- c(annotated_variants,
                              annotationEngine(variants, param, annotationCache,
                                               BPPARAM=BPPARAM))
    }

    message(sprintf("%d variants processed", n.var))
  }
  close(vcfFiles[[1]])

  gSO <- annotateSO(annotated_variants, sog(param))
  locMask <- conMask <- varTypMask <- logical(0)

  if (length(annotated_variants) > 0) {

    locMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$LOCATION)),
                                       levels(annotated_variants$LOCATION)))
    conMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$CONSEQUENCE)),
                                       levels(annotated_variants$CONSEQUENCE)))
    varTypMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$TYPE)),
                                          levels(annotated_variants$TYPE)))
  } else
    warning("No variants segregate following a de novo inheritance model.")

  MAFpopMask <- NA
  if (any(c("MafDb", "MafDb2") %in% sapply(param$otherAnnotations, class))) {
    ## assume AF columns are those containing AF[A-Z]+ and being of class 'numeric'
    cnAF <- colnames(mcols(annotated_variants))
    colsclasses <- sapply(mcols(annotated_variants), class)
    cnAF <- cnAF[intersect(grep("AF[A-Z]+", cnAF), grep("numeric", colsclasses))]
    MAFpopMask <- rep(TRUE, length(cnAF))
    names(MAFpopMask) <- cnAF
  }

  new("VariantFilteringResults", callObj=callobj, callStr=callstr, inputParameters=param,
      activeSamples=sampleNames, inheritanceModel="de novo", variants=annotated_variants,
      bamViews=BamViews(), gSO=gSO, filters=filters(param), cutoffs=cutoffs(param),
      dbSNPflag=NA_character_, OMIMflag=NA_character_,
      variantTypeMask=varTypMask, locationMask=locMask, consequenceMask=conMask, aaChangeType="Any",
      MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minScore5ss=NA_real_, minScore3ss=NA_real_, minCUFC=0)
})

.deNovoMask <- function(vObj, pedDf,
                        use=c("everything", "complete.obs", "all.obs"),
                        penetrance=1) {

  use <- match.arg(use)

  if (class(vObj) != "VRanges" && class(vObj) != "CollapsedVCF")
    stop("Argument 'vObj' should be either a 'VRanges' or a 'CollapsedVCF' object.")

  stopifnot(all(colnames(pedDf) %in% c("FamilyID", "IndividualID", "FatherID", "MotherID", "Gender", "Phenotype"))) ## QC

  nsamples <- nvariants <- 0
  if (class(vObj) == "VRanges") {
    nsamples <- nlevels(sampleNames(vObj))
    nvariants <- length(vObj)
  } else if (class(vObj) == "CollapsedVCF") {
    nsamples <- as.integer(ncol(vObj))
    nvariants <- nrow(vObj)
  }

  ## PENETRANCE ??

  if (sum(pedDf$Phenotype == 1) != 2 || sort(pedDf$Gender[pedDf$Phenotype == 1]) != 1:2) 
    stop("Current 'de novo' analysis requires two unaffected parents and one or more affected children.")

  unaff <- pedDf[pedDf$Phenotype == 1, ]
  unaff_dad <- unaff[unaff$Gender == 1, "IndividualID"]
  unaff_mom <- unaff[unaff$Gender == 2, "IndividualID"]
  aff <- pedDf[pedDf$Phenotype == 2, ]
  
  if (aff$FatherID != unaff_dad || aff$MotherID != unaff_mom)
    stop("Current 'de novo' analysis requires two unaffected parents and one or more affected children.")

  ## fetch genotypes
  gt <- NULL
  if (class(vObj) == "VRanges")
    gt <- do.call("cbind", split(vObj$GT, sampleNames(vObj)))
  else if (class(vObj) == "CollapsedVCF")
    gt <- geno(vObj)$GT

  missingMask <- apply(gt, 1, function(x) any(x == "." | x == "./." | x == ".|."))

  if (any(missingMask) && use == "all.obs")
    stop("There are missing genotypes and current policy to deal with them is 'all.obs', which does not allow them.")

  unaffectedMask <- rep(TRUE, times=nrow(gt))
  unaffgt <- gt[, unaff$IndividualID, drop=FALSE]
  if (any(missingMask) && use == "everything")
    unaffgt[unaffgt == "." | unaffgt == "./." | unaffgt == ".|."] <- NA_character_
  unaffectedMask <- unaffgt == "0/0" | unaffgt == "0|0"
  unaffectedMask <- apply(unaffectedMask, 1, all)
  rm(unaffgt)

  affgt <- gt[, aff$IndividualID, drop=FALSE]
  if (any(missingMask) && use == "everything")
    affgt[affgt == "." | affgt == "./." | affgt == ".|."] <- NA_character_
  affectedMask <- affgt == "0/1" | affgt == "0|1" | affgt == "1/1" | affgt == "1|1"
  affectedMask <- rowSums(affectedMask) == nrow(aff)
  rm(affgt)

  uaMask <- unaffectedMask & affectedMask
  if (any(missingMask) && use == "complete.obs")
    uaMask <- uaMask & !missingMask

  ## variants ultimately set to NA are discarded (should this be tuned by an argument?)
  uaMask[is.na(uaMask)] <- FALSE

  ## build logical mask for variants that segregate as a de novo trait
  denovoMask <- vector(mode="logical", length=nvariants) ## assume default values are FALSE

  if (class(vObj) == "VRanges") {
    idx <- split(1:nvariants, sampleNames(vObj))
    denovoMask[unlist(idx, use.names=FALSE)] <- rep(uaMask, times=nsamples)
  } else if (class(vObj) == "CollapsedVCF")
    denovoMask <- uaMask
  else
    warning(paste(sprintf("object 'vObj' has class %s, unknown to this function.",
                          class(vObj)),
                  "As a consequence, no variants are selected as do novo."))

  denovoMask
}

.deNovoFilter <- function(x) {

  if (is.null(param(x)$pedFilename))
    stop("Please specify a PED file name in the 'VariantFiltering' parameter object.")

  pedDf <- .readPEDfile(param(x)$pedFilename)

  .deNovoMask(vObj=allVariants(x, groupBy="nothing"),
              pedDf=pedDf, use="everything")
}
