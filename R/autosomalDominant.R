setMethod("autosomalDominant", signature(param="VariantFilteringParam"),
          function(param, BPPARAM=bpparam()) {

  ## store call for reproducing it later
  callobj <- match.call()
  callstr <- deparse(callobj)

  ## fetch necessary parameters
  vcfFiles <- param$vcfFiles
  ped <- param$pedFilename
  seqInfos <- param$seqInfos
  txdb <- param$txdb
  bsgeonme <- param$bsgenome

  if (!exists(as.character(substitute(BPPARAM))))
    stop(sprintf("Parallel back-end function %s given in argument 'BPPARAM' does not exist in the current workspace. Either you did not write correctly the function name or you did not load the package 'BiocParallel'.", as.character(substitute(BPPARAM))))

  if (length(vcfFiles) > 1)
    stop("More than one input VCF file is currently not supported. Please either merge the VCF files into a single one with software such as vcftools or GATK, or do the variant calling simultaneously on all samples, or proceed analyzing each file separately.")
  else if (length(vcfFiles) < 1)
    stop("A minimum of 1 vcf file has to be provided")

  pedf <- read.table(ped, header=FALSE, stringsAsFactors=FALSE)
  pedf <- pedf[, 1:6]
  colnames(pedf) <- c("FamilyID", "IndividualID", "FatherID", "MotherID", "Gender", "Phenotype")

  ## assuming Phenotype == 2 means affected and Phenotype == 1 means unaffected
  if (sum(pedf$Phenotype  == 2) < 1)
    stop("No affected individuals detected. Something is wrong with the PED file.")
  
  unaff <- pedf[ped$Phenotype == 1, ]
  aff <- pedf[pedf$Phenotype == 2, ]

  annotated_variants <- VRanges()
  open(vcfFiles[[1]])
  n.var <- 0
  while(nrow(vcf <- readVcf(vcfFiles[[1]], genome=seqInfos[[1]]))) {
  
    ## insert an index for each variant in the VCF file
    info(header(vcf)) <- rbind(info(header(vcf)),
                               DataFrame(Number=1, Type="Integer",
                                         Description="Variant index in the VCF file.",
                                         row.names="VCFIDX"))
    info(vcf)$VCFIDX <- (n.var+1):(n.var+nrow(vcf))
    varIDs <- names(rowData(vcf))

    n.var <- n.var + nrow(vcf)

    ## build logical masks of affected and unaffected individuals
    ## variants in unaffected individuals should be homozygous reference and
    ## in affected individuals should be either homozygous alternative or heterozygous alternative
    unaffectedMask <- rep(TRUE, times=nrow(vcf))
    if (nrow(unaff) > 0) {
      unaffectedMask <- geno(vcf)$GT[, unaff$IndividualID, drop=FALSE] == "0/0"
      unaffectedMask <- apply(unaffectedMask, 1, all)
    }

    affectedMask <- geno(vcf)$GT[, aff$IndividualID, drop=FALSE] == "0/1" |
                    geno(vcf)$GT[, aff$IndividualID, drop=FALSE] == "1/1"

    ## filter out variants that do not segregate as a "de novo" trait
    vcf <- vcf[unaffectedMask & affectedMask, ]

    ## coerce the VCF object to a VRanges object
    variants <- as(vcf, "VRanges")

    ## since the conversion of VCF to VRanges strips the VCF ID field, let's put it back
    variants$VARID <- varIDs[variants$VCFIDX]

    ## harmonize Seqinfo data between variants, annotations and reference genome
    variants <- .matchSeqinfo(variants, txdb, bsgenome)
  
    ## annotate variants
    annotated_variants <- c(annotated_variants, annotationEngine(variants, param, BPPARAM=BPPARAM))
  }
  close(vcfFiles[[1]])

  locMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$LOCATION)),
                                     levels(annotated_variants$LOCATION)))
  conMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$CONSEQUENCE)),
                                     levels(annotated_variants$CONSEQUENCE)))
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
      inheritanceModel="autosomal dominant", variants=dominant_annotated,
      dbSNPflag=NA_character_, OMIMflag=NA_character_, variantType="Any",
      locationMask=locMask, consequenceMask=conMask, aaChangeType="Any",
      MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minCRYP5ss=NA_real_, minCRYP3ss=NA_real_, minCUFC=0)
})
