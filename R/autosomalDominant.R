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
  txdb <- get(param$txdb)
  bsgenome <- get(param$bsgenome)
  sampleNames <- param$sampleNames

  if (!exists(as.character(substitute(BPPARAM))[1]))
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
  metadata(mcols(annotated_variants)) <- list(cutoffs=CutoffsList(), sortings=CutoffsList())

  open(vcfFiles[[1]])
  n.var <- 0
  flag <- TRUE
  while(flag && nrow(vcf <- readVcf(vcfFiles[[1]], genome=seqInfos[[1]], param=svparam))) {
  
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
      if (length(annotated_variants) > 0)
        annotated_variants <- c(annotated_variants,
                                annotationEngine(variants, param, annotationCache,
                                                 BPPARAM=BPPARAM))
      else
        annotated_variants <- annotationEngine(variants, param, annotationCache,
                                               BPPARAM=BPPARAM)

      if (length(vcfWhich(svparam)) > 0) ## temporary fix to keep this looping
        flag <- FALSE                    ## structure with access through genomic ranges
    }

    message(sprintf("%d variants processed", n.var))
  }
  close(vcfFiles[[1]])

  gSO <- annotateSO(annotated_variants, sog(param))
  annotated_variants <- addSOmetadata(annotated_variants)

  if (length(annotated_variants) == 0)
    warning("No variants segregate following an autosomal dominant inheritance model.")

  annoGroups <- list()
  if (!is.null(mcols(mcols(annotated_variants))$TAB)) {
    mask <- !is.na(mcols(mcols(annotated_variants))$TAB)
    annoGroups <- split(colnames(mcols(annotated_variants))[mask],
                      mcols(mcols(annotated_variants))$TAB[mask])
  }

  ## add functional annotation filters generated by the annotation engine
  funFilters <- FilterRules(lapply(metadata(mcols(annotated_variants))$filters,
                                   function(f) { environment(f) <- baseenv() ; f}))
  active(funFilters) <- FALSE ## by default, functional annotation filters are inactive
  flt <- c(filters(param), funFilters)
  fltMd <- rbind(filtersMetadata(param),
                 DataFrame(Description=sapply(metadata(mcols(annotated_variants))$filters, 
                                              attr, "description"),
                           AnnoGroup=sapply(metadata(mcols(annotated_variants))$filters,
                                            attr, "TAB")))
  cutoffs <- metadata(mcols(annotated_variants))$cutoffs
  sortings <- metadata(mcols(annotated_variants))$sortings
  bsgenome <- get(param$bsgenome)

  new("VariantFilteringResults", callObj=callobj, callStr=callstr,
      genomeDescription=as(bsgenome, "GenomeDescription"), inputParameters=param,
      activeSamples=sampleNames, inheritanceModel="autosomal dominant",
      variants=annotated_variants, bamViews=BamViews(), gSO=gSO, filters=flt,
      filtersMetadata=fltMd, cutoffs=cutoffs, sortings=sortings, annoGroups=annoGroups,
      minScore5ss=NA_real_, minScore3ss=NA_real_)
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

  stopifnot(all(colnames(pedDf) %in% c("FamilyID", "IndividualID", "FatherID", "MotherID", "Sex", "Phenotype"))) ## QC

  nsamples <- nvariants <- 0
  if (class(vObj) == "VRanges") {
    nsamples <- nlevels(sampleNames(vObj))
    nvariants <- length(vObj)
  } else if (class(vObj) == "CollapsedVCF") {
    nsamples <- as.integer(ncol(vObj))
    nvariants <- nrow(vObj)
  }

  ## assuming Phenotype == 2 means affected and Phenotype == 1 means unaffected
  if (sum(pedDf$Phenotype  == 2) < 1)
    stop("No affected individuals found in the PED file.")

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

  ## further restrict affected and unaffected individuals to
  ## those who have been genotyped
  gtind <- colnames(gt)
  unaff <- unaff[unaff$IndividualID %in% gtind, , drop=FALSE]
  aff <- aff[aff$IndividualID %in% gtind, , drop=FALSE]
  if (nrow(aff) == 0)
    stop("No affected individuals have genotypes.")

  phasedgt <- any(grepl("\\|", gt))

  missingMask <- rowSums(gt == ".") > 0
  if (phasedgt)
    missingMask <- missingMask | (rowSums(gt == ".|.") > 0)
  else
    missingMask <- missingMask | (rowSums(gt == "./.") > 0)

  ## missingMask <- apply(gt, 1, function(x) any(x == "." | x == "./." | x == ".|."))

  if (any(missingMask) && use == "all.obs")
    stop("There are missing genotypes and current policy to deal with them is 'all.obs', which does not allow them.")

  ## build logical masks of unaffected and affected individuals
  ## unaffected individuals should be homozygous reference and affected
  ## should be heterozygous of homozygous alternative
  unaffectedMask <- rep(TRUE, times=nrow(gt))
  if (nrow(unaff) > 0) {
    unaffgt <- gt[, unaff$IndividualID, drop=FALSE]
    if (any(missingMask) && use == "everything") {
      unaffgt[unaffgt == "."] <- NA_character_
      if (phasedgt)
        unaffgt[unaffgt == ".|."] <- NA_character_
      else
        unaffgt[unaffgt == "./."] <- NA_character_
    }
    if (phasedgt)
      unaffectedMask <- unaffgt == "0|0"
    else
      unaffectedMask <- unaffgt == "0/0"
    unaffectedMask <- rowSums(unaffectedMask, na.rm=TRUE) == nrow(unaff)
    rm(unaffgt)
  }

  affgt <- gt[, aff$IndividualID, drop=FALSE]
  if (any(missingMask) && use == "everything") {
    affgt[affgt == "."] <- NA_character_
    if (phasedgt)
      affgt[affgt == ".|."] <- NA_character_
    else
      affgt[affgt == "./."] <- NA_character_
  }

  ## old mask assuming biallelic variants
  ## affectedMask <- affgt == "0/1" | affgt == "0|1" | affgt == "1/1" | affgt == "1|1"
  ## affected mask is built taking into account multiallelic variants that typically happen in indels,
  ## i.e., 0/1, 0/2, 0/3, ..., 1/1, 2/2, 3/3, etc.
  if (phasedgt) {
    affectedMask <- matrix(grepl("0\\|\\|0", affgt), nrow=nrow(affgt)) & affgt != "0|0"
    affgt[affgt == "0|0"] <- NA_character_
    affgt <- strsplit(affgt, "|")
  } else {
    affectedMask <- matrix(grepl("0/", affgt), nrow=nrow(affgt)) & affgt != "0/0"
    affgt[affgt == "0/0"] <- NA_character_
    affgt <- strsplit(affgt, "/")
  }

  naMask <- sapply(affgt, function(x) any(is.na(x)))
  affectedMask <- as.vector(affectedMask)
  affectedMask[!naMask] <- affectedMask[!naMask] |
                           as.vector(tapply(unlist(affgt[!naMask]), rep(1:sum(!naMask), each=2),
                                            function(x) all(x[1] == x[2])))
  affectedMask <- matrix(affectedMask, ncol=nrow(aff))

  affectedMask <- rowSums(affectedMask, na.rm=TRUE) == nrow(aff)
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
