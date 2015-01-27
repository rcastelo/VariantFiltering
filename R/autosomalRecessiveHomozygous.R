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
  
  pedf <- read.table(ped, header=FALSE)

  if (sum(pedf[, 6] == 2) < 1)
    stop("No affected individuals detected. Something is wrong with the PED file.")
  
  unaff <- pedf[pedf[, 6] == 1, ]
  aff <- pedf[pedf[, 6] == 2, ]
  
  annotated_variants <- GRanges()
  open(vcfFiles[[1]])
  while (nrow(vcf <- readVcf(vcfFiles[[1]], genome=seqInfos[[1]]))) {
    message(sprintf("%d variants read from input VCF file into main memory.", nrow(vcf)))
   
    carriers <- switch(nrow(unaff),
                       one_ind_ms(vcf, "0/1", unaff, filterTag),
                       two_ind_ms(vcf, "0/1", unaff, filterTag),
                       three_ind_ms(vcf, "0/1", unaff, filterTag),
                       four_ind_ms(vcf, "0/1", unaff, filterTag),
                       five_ind_ms(vcf, "0/1", unaff, filterTag))
  
    affected <- switch(nrow(aff)+1,
                       stop("No affected individuals detected. Something is wrong with the PED file."),
                       one_ind_ms(vcf, "1/1", aff, filterTag),
                       two_ind_ms(vcf, "1/1", aff, filterTag),
                       three_ind_ms(vcf, "1/1", aff, filterTag),
                       four_ind_ms(vcf, "1/1", aff, filterTag),
                       five_ind_ms(vcf, "1/1", aff, filterTag))
  
    variants <- affected
    if (length(carriers) >= 1) {
      variants <- affected[sharedVariants(affected, carriers)]
    }
  
    variants <- .matchSeqinfo(variants, txdb, bsgenome)

    annotated_variants <- c(annotated_variants, annotationEngine(variants, param, BPPARAM=BPPARAM))
  }
  close(vcfFiles[[1]])

  MAFpopMask <- NA
  if ("MafDb" %in% sapply(param$otherAnnotations, class)) {
    cnAF <- colnames(mcols(annotated_variants))
    cnAF <- cnAF[grep("AF", cnAF)]
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
