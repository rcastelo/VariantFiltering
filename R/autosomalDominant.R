setMethod("autosomalDominant", signature(param="VariantFilteringParam"),
          function(param, svparam=ScanVcfParam(),
                   use=c("everything", "complete.obs", "all.obs"),
                   BPPARAM=bpparam("SerialParam")) {

  ## store call for reproducing it later
  callobj <- match.call()
  callstr <- gsub(".local", "autosomalDominant", deparse(callobj))

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

  unaff <- pedDf[pedDf$Phenotype == 1, ]
  aff <- pedDf[pedDf$Phenotype == 2, ]

  annotationCache <- new.env() ## cache annotations when using VariantAnnotation::locateVariants()
  annotated_variants <- VRanges()
  open(vcfFiles[[1]])
  n.var <- 0
  while(nrow(vcf <- readVcf(vcfFiles[[1]], genome=seqInfos[[1]], param=svparam))) {
  
    ## insert an index for each variant in the VCF file
    info(header(vcf)) <- rbind(info(header(vcf)),
                               DataFrame(Number=1, Type="Integer",
                                         Description="Variant index in the VCF file.",
                                         row.names="VCFIDX"))
    info(vcf)$VCFIDX <- (n.var+1):(n.var+nrow(vcf))
    varIDs <- rownames(vcf)

    n.var <- n.var + nrow(vcf)

    mask <- .autosomalDominantMask(vcf, pedDf, bsgenome, use)

    if (any(mask)) {

      ## filter out variants that do not segregate as an autosomal dominant trait
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
    warning("No variants segregate following an autosomal dominant inheritance model.")

  MAFpopMask <- NA
  if ("MafDb" %in% sapply(param$otherAnnotations, class)) {
    ## assume AF columns are those containing AF[A-Z]+ and being of class 'numeric'
    cnAF <- colnames(mcols(annotated_variants))
    colsclasses <- sapply(mcols(annotated_variants), class)
    cnAF <- cnAF[intersect(grep("AF[A-Z]+", cnAF), grep("numeric", colsclasses))]
    MAFpopMask <- rep(TRUE, length(cnAF))
    names(MAFpopMask) <- cnAF
  }

  new("VariantFilteringResults", callObj=callobj, callStr=callstr, inputParameters=param,
      activeSamples=sampleNames, inheritanceModel="autosomal dominant",
      variants=annotated_variants, bamViews=BamViews(), gSO=gSO, filters=filters(param),
      cutoffs=cutoffs(param), dbSNPflag=NA_character_, OMIMflag=NA_character_,
      variantTypeMask=varTypMask, locationMask=locMask, consequenceMask=conMask, aaChangeType="Any",
      MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minScore5ss=NA_real_, minScore3ss=NA_real_, minCUFC=0)
})

## build a logical mask whose truth values correspond to variants that segregate
## according to an autosomal dominant inheritance model: variants in unaffected
## individuals should be homozygous reference and in affected individuals should
## be either homozygous alternative or heterozygous
.autosomalDominantMask <- function(vObj, pedDf, bsgenome,
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

  ## assuming Phenotype == 2 means affected and Phenotype == 1 means unaffected
  if (sum(pedDf$Phenotype  == 2) < 1)
    stop("No affected individuals detected. Something is wrong with the PED file.")

  unaff <- pedDf[pedDf$Phenotype == 1, ]
  aff <- pedDf[pedDf$Phenotype == 2, ]

  ## restrict upfront variants to those in autosomal chromosomes
  ## we subset to the first element of the value returned by seqlevelsStyle()
  ## to deal with cases in which only a subset of chromosomes is contained in
  ## the input VCF (typically for teaching/example/illustration purposes) which
  ## matches more than one chromosome style, or because Ensembl is identical to NCBI for human :\
  snames <- as.character(seqnames(vObj))
  autosomalMask <- snames %in% extractSeqlevelsByGroup(organism(bsgenome),
                                                       seqlevelsStyle(vObj)[1],
                                                       group="auto")

  ## build logical mask for variants that segregate as an autosomal dominant trait
  adomMask <- vector(mode="logical", length=nvariants) ## assume default values are FALSE

  if (!any(autosomalMask))
    return(adomMask)

  ## fetch genotypes
  gt <- NULL
  if (class(vObj) == "VRanges")
    gt <- do.call("cbind", split(vObj$GT[autosomalMask], sampleNames(vObj)))
  else if (class(vObj) == "CollapsedVCF")
    gt <- geno(vObj)$GT[autosomalMask, , drop=FALSE]

  missingMask <- apply(gt, 1, function(x) any(x == "." | x == "./." | x == ".|."))

  if (any(missingMask) && use == "all.obs")
    stop("There are missing genotypes and current policy to deal with them is 'all.obs', which does not allow them.")

  unaffectedMask <- rep(TRUE, times=nrow(gt))
  if (nrow(unaff) > 0) {
    unaffgt <- gt[, unaff$IndividualID, drop=FALSE]
    if (any(missingMask) && use == "everything")
      unaffgt[unaffgt == "." | unaffgt == "./." | unaffgt == ".|."] <- NA_character_
    unaffectedMask <- unaffgt == "0/0" | unaffgt == "0|0"
    unaffectedMask <- apply(unaffectedMask, 1, all)
    rm(unaffgt)
  }

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

  if (class(vObj) == "VRanges") {
    nauto <- sum(autosomalMask)
    idx <- split(1:nauto, sampleNames(vObj[autosomalMask]))
    mask <- vector(mode="logical", length=nauto)
    mask[unlist(idx, use.names=FALSE)] <- rep(uaMask, times=nsamples)
    adomMask[autosomalMask] <- mask
  } else if (class(vObj) == "CollapsedVCF")
    adomMask[autosomalMask] <- uaMask
  else
    warning(paste(sprintf("object 'vObj' has class %s, unknown to this function.",
                          class(vObj)),
                  "As a consequence, no variants are selected as autosomal dominant."))

  adomMask
}

.autosomalDominantFilter <- function(x) {

  if (is.null(param(x)$pedFilename))
    stop("Please specify a PED file name in the 'VariantFiltering' parameter object.")

  pedDf <- .readPEDfile(param(x)$pedFilename)

  .autosomalDominantMask(vObj=allVariants(x, groupBy="nothing"), pedDf=pedDf,
                         bsgenome=param(x)$bsgenome, use="everything")
}
