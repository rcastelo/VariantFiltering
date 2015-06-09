setMethod("autosomalRecessiveHomozygous", signature(param="VariantFilteringParam"),
          function(param, BPPARAM=bpparam("SerialParam")) {

  ## store call for reproducing it later
  callobj <- match.call()
  callstr <- gsub(".local", "autosomalRecessiveHomozygous", deparse(callobj))

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
    stop("A minimum of 1 VCF file has to be provided.")
  
  if (is.na(ped))
    stop("Please specify a PED file name when building the parameter object.")

  if (!file.exists(ped))
    stop(sprintf("could not open the PED file %s.", ped))

  pedf <- read.table(ped, header=FALSE, stringsAsFactors=FALSE)
  pedf <- pedf[, 1:6]
  colnames(pedf) <- c("FamilyID", "IndividualID", "FatherID", "MotherID", "Gender", "Phenotype")

  ## assuming Phenotype == 2 means affected and Phenotype == 1 means unaffected
  if (sum(pedf$Phenotype  == 2) < 1)
    stop("No affected individuals detected. Something is wrong with the PED file.")
  
  unaff <- pedf[pedf$Phenotype == 1, ]
  aff <- pedf[pedf$Phenotype == 2, ]
  
  annotationCache <- new.env() ## cache annotations when using VariantAnnotation::locateVariants()
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
    varIDs <- rownames(vcf)

    n.var <- n.var + nrow(vcf)

    ## restrict upfront variants to those in autosomal chromosomes
    autosomalMask <- seqnames(vcf) %in% extractSeqlevelsByGroup(organism(bsgenome),
                                                                seqlevelsStyle(vcf),
                                                                group="auto")
    vcf <- vcf[autosomalMask, ]

    ## build logical masks of carriers (unaffected) and affected individuals
    ## variants in carriers should be heterozygous and affected should be homozygous alternative
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
    annotated_variants <- c(annotated_variants, annotationEngine(variants, param, annotationCache,
                                                                 BPPARAM=BPPARAM))

    message(sprintf("%d variants processed", n.var))
  }
  close(vcfFiles[[1]])

  gSO <- annotateSO(annotated_variants, sog(param))

  locMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$LOCATION)),
                                     levels(annotated_variants$LOCATION)))
  conMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$CONSEQUENCE)),
                                     levels(annotated_variants$CONSEQUENCE)))
  varTypMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$TYPE)),
                                        levels(annotated_variants$TYPE)))
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
      activeSamples=sampleNames, inheritanceModel="autosomal recessive homozygous",
      variants=annotated_variants, bamViews=BamViews(), gSO=gSO, filters=filters(param),
      cutoffs=cutoffs(param), dbSNPflag=NA_character_, OMIMflag=NA_character_,
      variantTypeMask=varTypMask, locationMask=locMask, consequenceMask=conMask, aaChangeType="Any",
      MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minScore5ss=NA_real_, minScore3ss=NA_real_, minCUFC=0)
})
