setMethod("autosomalRecessiveHomozygous", signature(param="VariantFilteringParam"),
          function(param, BPPARAM=bpparam()) {

  callobj <- match.call()
  callstr <- deparse(callobj)
  vcfFiles <- param$vcfFiles
  ped <- param$pedFilename
  seqInfos <- param$seqInfos
  txdb <- param$txdb
  bsgenome <- param$bsgenome
  filterTag <- param$filterTag
    
  if (!exists(as.character(substitute(BPPARAM))))
    stop(sprintf("Parallel back-end function %s given in argument 'BPPARAM' does not exist in the current workspace. Either you did not write correctly the function name or you did not load the package 'BiocParallel'.", as.character(substitute(BPPARAM))))
  
  if (length(vcfFiles) > 1) {
    stop("More than one input VCF file is currently not supported. Please either merge the VCF files into a single one with vcftools, do the variant calling simultaneously on all samples, or proceed analysing each file separately.")
  } else if (length(vcfFiles) < 1)
    stop("A minimum of 1 vcf file has to be provided")
  
  pedf <- read.table(ped, header=FALSE, stringsAsFactors=FALSE)
  pedf <- pedf[, 1:6]
  colnames(pedf) <- c("FamilyID", "IndividualID", "FatherID", "MotherID", "Gender", "Phenotype")

  ## assuming Phenotype == 2 means affected and Phenotype == 1 means unaffected
  if (sum(pedf$Phenotype  == 2) < 1)
    stop("No affected individuals detected. Something is wrong with the PED file.")
  
  unaff <- pedf[pedf$Phenotype == 1, ]
  aff <- pedf[pedf$Phenotype == 2, ]
  
  annotated_variants <- VRanges()
  open(vcfFiles[[1]])
  n.var <- 0
  while (nrow(vcf <- readVcf(vcfFiles[[1]], genome=seqInfos[[1]]))) {
   
    ## insert an index for each variant in the VCF file
    info(header(vcf)) <- rbind(info(header(vcf)),
                               DataFrame(Number=1, Type="Integer",
                                         Description="Variant index in the VCF file.",
                                         row.names="VCFIDX"))
    info(vcf)$VCFIDX <- (n.var+1):(n.var+nrow(vcf))
    varIDs <- names(rowData(vcf))

    n.var <- n.var + nrow(vcf)

    carriersMask <- rep(TRUE, times=nrow(vcf))
    if (nrow(unaff) > 0) {
      carriersMask <- geno(vcf)$GT[, unaff$IndividualID, drop=FALSE] == "0/1"
      carriersMask <- apply(carriersMask, 1, all)
    }

    affectedMask <- geno(vcf)$GT[, aff$IndividualID, drop=FALSE] == "1/1"
    affectedMask <- apply(affectedMask, 1, all)

    ## filter out variants that do not segregate as an autosomal recessive homozygous trait
    vcf <- vcf[carriersMask & affectedMask, ]

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
    cnAF <- cnAF[intersect(grep("AF[A-Z]+", cnAF), grep("numeric", colsclasses))]
    MAFpopMask <- rep(TRUE, length(cnAF))
    names(MAFpopMask) <- cnAF
  }

  new("VariantFilteringResults", callObj=callobj, callStr=callstr, inputParameters=param,
      inheritanceModel="autosomal recessive homozygous", variants=annotated_variants,
      dbSNPflag=NA_character_, OMIMflag=NA_character_, variantType="Any",
      aaChangeType="Any", MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minCRYP5ss=NA_real_, minCRYP3ss=NA_real_)
})
