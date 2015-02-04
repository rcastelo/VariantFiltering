setMethod("unrelatedIndividuals", signature(param="VariantFilteringParam"),
          function(param, BPPARAM=bpparam()) {

  callobj <- match.call()
  callstr <- deparse(callobj)
  vcfFiles <- param$vcfFiles
  seqInfos <- param$seqInfos
  txdb <- param$txdb
  bsgenome <- param$bsgenome
  
  if (!exists(as.character(substitute(BPPARAM))))
    stop(sprintf("Parallel back-end function %s given in argument 'BPPARAM' does not exist in the current workspace. Either you did not write correctly the function name or you did not load the package 'BiocParallel'.", as.character(substitute(BPPARAM))))
  
  if (length(vcfFiles) > 1) {
    stop("More than one input VCF file is currently not supported. Please either merge the VCF files into a single one with vcftools, do the variant calling simultaneously on all samples, or proceed analyzing each file separately.")
  } else if (length(vcfFiles) < 1)
    stop("A minimum of 1 vcf file has to be provided")

  
  annotated_variants <- VRanges()
  open(vcfFiles[[1]])
  n.var <- 0
  while(nrow(vcf <- readVcf(vcfFiles[[1]], genome=seqInfos[[1]]))) {

    # insert an index for each variant in the VCF file
    info(header(vcf)) <- rbind(info(header(vcf)),
                               DataFrame(Number=1, Type="Integer",
                                         Description="Variant index in the VCF file.",
                                         row.names="VCFIDX"))
    info(vcf)$VCFIDX <- (n.var+1):(n.var+nrow(vcf))
    varIDs <- names(rowData(vcf))

    n.var <- n.var + nrow(vcf)

    ## coerce the VCF object to a VRanges object
    variants <- as(vcf, "VRanges")

    ## since the conversion of VCF to VRanges strips the VCF ID field, let's put it back
    variants$VARID <- varIDs[variants$VCFIDX]

    ## harmonize Seqinfo data between variants, annotations and reference genome
    variants <- .matchSeqinfo(variants, txdb, bsgenome)

    ## annotate variants
    annotated_variants <- c(annotated_variants, annotationEngine(variants, param, BPPARAM=BPPARAM))

    message(sprintf("%d variants processed", n.var))
  }
  close(vcfFiles[[1]])

  MAFpopMask <- NA
  if ("MafDb" %in% sapply(param$otherAnnotations, class)) {
    ## assume AF columns are those containing AF[A-Z]+ and being of class 'numeric'
    cnAF <- colnames(mcols(annotated_variants))
    colsclasses <- sapply(mcols(annotated_variants), class)
    nAF <- cnAF[intersect(grep("AF[A-Z]+", cnAF), grep("numeric", colsclasses))]
    MAFpopMask <- rep(TRUE, length(cnAF))
    names(MAFpopMask) <- cnAF
  }

  new("VariantFilteringResultsUI", callObj=callobj, callStr=callstr, inputParameters=param,
      inheritanceModel="unrelated individuals", variants=annotated_variants,
      indselected=NA_character_, selectgene=NA_character_, dbSNPflag=NA_character_, OMIMflag=NA_character_,
      variantType="Any", aaChangeType="Any", MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minCRYP5ss=NA_real_, minCRYP3ss=NA_real_)
})
